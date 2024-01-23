import os

OGS_USE_PATH = os.getenv("OGS_USE_PATH", "False").lower() in ("true", "1", "t")
