#!/usr/bin/env python3
"""Verify the frozen cap-13 certificate from the portable bundle."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHECKER = ROOT / "verification" / "verify_certificate.py"
CERTIFICATE = (
    ROOT / "artifacts" / "certificates" / "d15-cap13-certificate.json"
)


def main() -> int:
    for label, path in (("checker", CHECKER), ("certificate", CERTIFICATE)):
        if not path.is_file():
            raise SystemExit(f"{label} is missing from the bundle: {path.relative_to(ROOT)}")

    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONHASHSEED": "0",
            "PYTHONUNBUFFERED": "1",
        }
    )
    command = [
        sys.executable,
        str(CHECKER),
        str(CERTIFICATE),
        "--mutation-tests",
        "--verbose",
    ]
    return subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        check=False,
    ).returncode


if __name__ == "__main__":
    raise SystemExit(main())
