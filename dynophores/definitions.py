"""
Definitions for dynophore generation and visualization.
"""

from pathlib import Path

DYNOPHORE_JAR_NAME = "dynophore.jar"
DYNOPHORE_JAR_URL = f"https://drug-design.de/dynophore/downloads/{DYNOPHORE_JAR_NAME}"
DYNOPHORE_JAR_PATH = Path(__file__).parent / "generate" / DYNOPHORE_JAR_NAME
