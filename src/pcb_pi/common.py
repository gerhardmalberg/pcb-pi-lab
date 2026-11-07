"""Small shared utilities; no numerical dependencies on the host."""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

class WorkflowError(RuntimeError):
    pass

def digest(path):
    with open(path, 'rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)

def number(value, name, minimum=0, allow_zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise WorkflowError(f'{name} must be a finite number')
    if value < minimum or (value == minimum and not allow_zero):
        raise WorkflowError(f'{name} must be {">=" if allow_zero else ">"} {minimum}')
    return float(value)
