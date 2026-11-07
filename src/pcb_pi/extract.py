"""Read KiCad using its own Python bindings. Never save or refill the board."""
import argparse
import json
import math
import os
import re
from pathlib import Path
from .common import WorkflowError, digest, write_json


def expression(text, keyword):
    start = text.index("(" + keyword)
    depth = 0
    quoted = escaped = False
    for i in range(start, len(text)):
        c = text[i]
        if escaped:
            escaped = False
        elif c == "\\" and quoted:
            escaped = True
        elif c == '"':
            quoted = not quoted
        elif not quoted:
            depth += (c == "(") - (c == ")")
            if depth == 0:
                return text[start:i+1]
    raise WorkflowError(f"Unclosed KiCad expression: {keyword}")


def sexpr(text):
    tokens = re.findall(r'"(?:\\.|[^"\\])*"|[()]|[^\s()]+', text)
    stack = [[]]
    for token in tokens:
        if token == "(":
            stack.append([])
        elif token == ")":
            value = stack.pop()
            stack[-1].append(value)
        else:
            stack[-1].append(json.loads(token) if token.startswith('"') else token)
    return stack[0][0]


def stackup(text):
    if not re.search(r'\(stackup\s',text):
        raise WorkflowError('PCB has no saved physical stackup; configure layer thicknesses and materials in KiCad Board Setup, save the board, then rerun extraction')
    layers = []
    for entry in sexpr(expression(text, "stackup"))[1:]:
        if not isinstance(entry, list) or entry[0] != "layer":
            continue
        name = entry[1]
        pieces = [{}]
        for value in entry[2:]:
            if value == "addsublayer":
                pieces.append({"type": pieces[0].get("type", "dielectric")})
            elif isinstance(value, list) and len(value) >= 2:
                pieces[-1][value[0]] = value[1]
        for index, item in enumerate(pieces):
            thickness = float(item.get("thickness", 0))
            if thickness <= 0:
                continue
            layers.append({"name": name, "sublayer": index, "thickness_mm": thickness,
                           "kind": "copper" if name.endswith(".Cu") else "dielectric",
                           "material": item.get("material", item.get("type", name)),
                           "epsilon_r": float(item["epsilon_r"]) if "epsilon_r" in item else None,
                           "loss_tangent": float(item.get("loss_tangent", 0))})
    if not any(x["name"] == "F.Cu" for x in layers):
        raise WorkflowError("PCB has no physical F.Cu stackup; supply a saved stackup first")
    z = sum(x["thickness_mm"] for x in layers[:next(i for i,x in enumerate(layers) if x["name"] == "F.Cu")])
    for item in layers:
        item["z_top_mm"] = z
        z -= item["thickness_mm"]
        item["z_bottom_mm"] = z
    return layers


def polygons(poly):
    def ring(chain):
        return [[chain.CPoint(i).x / 1e6, chain.CPoint(i).y / 1e6] for i in range(chain.PointCount())]
    return [{"outer": ring(poly.COutline(i)),
             "holes": [ring(poly.CHole(i,j)) for j in range(poly.HoleCount(i))]}
            for i in range(poly.OutlineCount())]


def extract(board_path, output, chord_mm=0.002):
    import pcbnew
    board_path = Path(board_path).resolve()
    before = digest(board_path)
    board = pcbnew.LoadBoard(str(board_path))
    layers = stackup(board_path.read_text())
    copper_ids = {board.GetLayerID(x["name"]): x["name"] for x in layers if x["kind"] == "copper"}
    outline = pcbnew.SHAPE_POLY_SET()
    if not board.GetBoardPolygonOutlines(outline, False, None, False, False):
        raise WorkflowError("Board outline cannot be extracted as closed polygons")
    result = {"schema": 1, "source": board_path.name, "sha256": before,
              "kicad_version": pcbnew.GetBuildVersion(), "units": "mm",
              "coordinates": "saved KiCad XY; z=0 at F.Cu top, positive above board",
              "polygon_chord_error_mm": chord_mm, "stackup": layers,
              "declared_thickness_mm": board.GetDesignSettings().GetBoardThickness()/1e6,
              "outline": polygons(outline), "copper": [], "holes": [],
              "components": [], "keepouts": [], "warnings": []}
    # KiCad applies saved pad/footprint/board mask margins, tented vias and
    # mask-layer graphics here. Copper pad outlines alone are not mask openings.
    result['soldermask_openings'] = {}
    for name, lid in [('F.Mask', pcbnew.F_Mask), ('B.Mask', pcbnew.B_Mask)]:
        poly = pcbnew.SHAPE_POLY_SET()
        board.ConvertBrdLayerToPolygonalContours(lid, poly)
        result['soldermask_openings'][name] = polygons(poly)
    def copper(item, layer, kind, extra=None):
        # Layer span includes the plated barrel even when an unused annular
        # ring is removed. Ask KiCad which pad/via copper is actually flashed.
        if kind in ('pad', 'via') and not item.FlashLayer(layer):
            return
        poly = pcbnew.SHAPE_POLY_SET()
        item.TransformShapeToPolygon(poly, layer, 0, round(chord_mm*1e6), pcbnew.ERROR_INSIDE)
        if poly.OutlineCount():
            result["copper"].append({"kind": kind, "layer": copper_ids[layer],
                                     "net": item.GetNetname(), "polygons": polygons(poly), **(extra or {})})
    for track in board.GetTracks():
        if track.GetClass() == "PCB_VIA":
            top, bottom = track.TopLayer(), track.BottomLayer()
            enabled = [lid for lid in copper_ids if track.IsOnLayer(lid)]
            for lid in enabled:
                copper(track, lid, "via")
            result["holes"].append({"kind": "via", "center": [track.GetPosition().x/1e6,track.GetPosition().y/1e6],
                                    "drill": [track.GetDrillValue()/1e6]*2, "angle_deg": 0,
                                    "plated": True, "net": track.GetNetname(),
                                    "layers": [copper_ids[x] for x in enabled]})
        elif track.GetLayer() in copper_ids:
            copper(track, track.GetLayer(), "track", {"start": [track.GetStart().x/1e6,track.GetStart().y/1e6],
                   "end": [track.GetEnd().x/1e6,track.GetEnd().y/1e6], "width_mm": track.GetWidth()/1e6})
    for zone in board.Zones():
        if zone.GetIsRuleArea():
            result["keepouts"].append({"name":zone.GetZoneName(),
                                       "prohibits":{"tracks":zone.GetDoNotAllowTracks(),
                                           "vias":zone.GetDoNotAllowVias(),"pads":zone.GetDoNotAllowPads(),
                                           "zone_fills":zone.GetDoNotAllowZoneFills(),
                                           "footprints":zone.GetDoNotAllowFootprints()},
                                       "layers": [copper_ids[x] for x in copper_ids if zone.IsOnLayer(x)],
                                       "polygons": polygons(zone.Outline())})
            continue
        for lid in copper_ids:
            if zone.IsOnLayer(lid):
                if zone.HasFilledPolysForLayer(lid):
                    result["copper"].append({"kind": "zone", "layer": copper_ids[lid], "net": zone.GetNetname(),
                                             "polygons": polygons(zone.GetFilledPolysList(lid))})
                else:
                    result["warnings"].append(f"No saved fill for {zone.GetNetname()} on {copper_ids[lid]}")
    for footprint in board.GetFootprints():
        ref = footprint.GetReference()
        component = {"ref": ref, "value": footprint.GetValue(),
                     "footprint": str(footprint.GetFPID().GetLibItemName()),
                     "dnp": bool(footprint.IsDNP()),
                     "position": [footprint.GetPosition().x/1e6, footprint.GetPosition().y/1e6],
                     "angle_deg": footprint.GetOrientationDegrees(), "flipped": footprint.IsFlipped(),
                     "pads": [], "models": []}
        for pad in footprint.Pads():
            position = [pad.GetPosition().x/1e6,pad.GetPosition().y/1e6]
            enabled = [lid for lid in copper_ids if pad.IsOnLayer(lid)]
            component["pads"].append({"number": pad.GetNumber(), "net": pad.GetNetname(), "position": position,
                                      "size": [pad.GetSize().x/1e6,pad.GetSize().y/1e6],
                                      "angle_deg": pad.GetOrientationDegrees(),
                                      "layers": [copper_ids[x] for x in enabled]})
            if pad.GetAttribute() != pcbnew.PAD_ATTRIB_NPTH:
                for lid in enabled:
                    copper(pad, lid, "pad", {"ref": ref, "number": pad.GetNumber()})
            drill = [pad.GetDrillSize().x/1e6,pad.GetDrillSize().y/1e6]
            if max(drill) > 0:
                result["holes"].append({"kind": "pad", "center": position, "drill": drill,
                    "angle_deg": pad.GetOrientationDegrees(), "plated": pad.GetAttribute() != pcbnew.PAD_ATTRIB_NPTH,
                    "ref": ref, "number": pad.GetNumber(), "net": pad.GetNetname(),
                    "layers": [copper_ids[x] for x in enabled]})
        for model in footprint.Models():
            raw = model.m_Filename
            resolved = Path(os.path.expandvars(raw.replace("${KIPRJMOD}", str(board_path.parent))))
            if not resolved.is_absolute():resolved=board_path.parent/resolved
            resolved=resolved.resolve()
            component["models"].append({"path": raw, "resolved_path": str(resolved),
                "sha256": digest(resolved) if resolved.is_file() else None,
                "scale": [model.m_Scale.x, model.m_Scale.y, model.m_Scale.z],
                "offset": [model.m_Offset.x, model.m_Offset.y, model.m_Offset.z],
                "rotation": [model.m_Rotation.x, model.m_Rotation.y, model.m_Rotation.z]})
        result["components"].append(component)
    if digest(board_path) != before:
        raise WorkflowError("Board changed while extracting; rerun against a stable saved revision")
    total=sum(layer['thickness_mm'] for layer in layers)
    if abs(total-result['declared_thickness_mm'])>0.001:
        result['warnings'].append(f"Physical stackup sum {total:g} mm differs from board setting {result['declared_thickness_mm']:g} mm; meshing uses the physical layers")
    write_json(output, result)
    print(f"Extracted {len(result['copper'])} copper items, {len(result['components'])} components; source unchanged")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("board")
    p.add_argument("output")
    p.add_argument("--chord-mm", type=float, default=0.002)
    a = p.parse_args()
    extract(a.board,a.output,a.chord_mm)
