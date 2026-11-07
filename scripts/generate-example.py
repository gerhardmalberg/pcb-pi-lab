#!/usr/bin/env python3
"""Generate this repository's synthetic board only; no installed libraries needed."""
from pathlib import Path
import pcbnew as k

def present(board):
    """Center the fixture on landscape A4 and separate its reference labels."""
    for footprint in board.GetFootprints():
        label=footprint.Reference()
        label.SetLayer(k.F_SilkS)
        label.SetVisible(True)
        label.SetTextSize(k.VECTOR2I(k.FromMM(.9),k.FromMM(.9)))
        label.SetTextThickness(k.FromMM(.15))
        offset=k.VECTOR2I(k.FromMM(-1.5 if footprint.GetReference()=='U2' else 0),0)
        label.SetPosition(footprint.GetPosition()+offset)
        footprint.Value().SetVisible(False)
    shift=k.VECTOR2I(k.FromMM(148.5),k.FromMM(105))-board.GetBoardEdgesBoundingBox().GetCenter()
    for item in [*board.GetFootprints(),*board.GetTracks(),*board.GetDrawings(),*board.Zones()]:
        item.Move(shift)
    settings=board.GetDesignSettings()
    settings.SetAuxOrigin(settings.GetAuxOrigin()+shift)
    settings.SetGridOrigin(settings.GetGridOrigin()+shift)

def create(path):
    b=k.BOARD();b.SetCopperLayerCount(2);b.GetDesignSettings().SetBoardThickness(k.FromMM(1.036))
    p=lambda x,y:k.VECTOR2I(k.FromMM(x),k.FromMM(y))
    nets={}
    for name in ('SUPPLY','RETURN'):
        n=k.NETINFO_ITEM(b,name);b.Add(n);nets[name]=n
    corners=[(0,0),(20,0),(20,12),(0,12)]
    for a,z in zip(corners,corners[1:]+corners[:1]):
        e=k.PCB_SHAPE();e.SetShape(k.SHAPE_T_SEGMENT);e.SetStart(p(*a));e.SetEnd(p(*z));e.SetLayer(k.Edge_Cuts);e.SetWidth(k.FromMM(.05));b.Add(e)
    def track(a,z,net,width=.8,layer=k.F_Cu):
        t=k.PCB_TRACK(b);t.SetStart(p(*a));t.SetEnd(p(*z));t.SetWidth(k.FromMM(width));t.SetLayer(layer);t.SetNet(nets[net]);b.Add(t)
    def via(x,y,net):
        v=k.PCB_VIA(b);v.SetPosition(p(x,y));v.SetWidth(k.FromMM(.8));v.SetDrill(k.FromMM(.4));v.SetLayerPair(k.F_Cu,k.B_Cu);v.SetViaType(k.VIATYPE_THROUGH);v.SetNet(nets[net]);b.Add(v)
    for ref,x,y1,y2,value in [('J1',2,3,9,'Source'),('U1',18,3,9,'Load 1'),('U2',13,3,9,'Load 2'),('C1',16,4.5,7.5,'10u'),('C2',14,4.5,7.5,'1u')]:
        f=k.FOOTPRINT(b);f.SetReference(ref);f.SetValue(value);f.SetPosition(p(x,6))
        for num,y,net in [('1',y1,'SUPPLY'),('2',y2,'RETURN')]:
            pad=k.PAD(f);pad.SetNumber(num);pad.SetAttribute(k.PAD_ATTRIB_SMD);pad.SetShape(k.PAD_SHAPE_RECT);pad.SetSize(p(1,1));pad.SetPosition(p(x,y));pad.SetNet(nets[net]);ls=k.LSET();ls.AddLayer(k.F_Cu);pad.SetLayerSet(ls);f.Add(pad)
        b.Add(f)
        if ref.startswith('C'):
            track((x,3),(x,y1),'SUPPLY');track((x,y2),(x,9),'RETURN')
        via(x,9,'RETURN')
    track((2,3),(10,3),'SUPPLY',.25)
    via(10,3,'SUPPLY');via(11,3,'SUPPLY');track((10,3),(11,3),'SUPPLY',.8,k.B_Cu)
    track((11,3),(18,3),'SUPPLY')
    zone=k.ZONE(b);zone.SetLayer(k.B_Cu);zone.SetNet(nets['RETURN'])
    poly=k.SHAPE_POLY_SET();poly.NewOutline();zone.Outline().NewOutline()
    for xy in [(1,8),(19,8),(19,10),(1,10)]:poly.Append(p(*xy));zone.Outline().Append(p(*xy))
    zone.SetFilledPolysList(k.B_Cu,poly);zone.SetIsFilled(True);b.Add(zone)
    path=Path(path)
    if path.exists():raise RuntimeError('Refusing to overwrite an existing board')
    present(b)
    k.SaveBoard(str(path),b)
    stack='''(stackup
      (layer "F.Cu" (type "copper") (thickness 0.018))
      (layer "dielectric 1" (type "core") (thickness 1.0) (material "FR4") (epsilon_r 4.2))
      (layer "B.Cu" (type "copper") (thickness 0.018)))'''
    path.write_text(path.read_text().replace('(setup','(setup\n'+stack,1))

if __name__=='__main__':
    import sys
    create(sys.argv[1])
