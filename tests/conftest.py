"""Shared pytest setup: make the pipeline's scripts importable."""
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "app" / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
