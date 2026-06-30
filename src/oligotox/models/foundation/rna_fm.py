"""
RNA-FM integration stub for OligoTox.
Status: scaffolded automatically; replace with full implementation.
"""

class RNAFM:
    def __init__(self, config=None):
        self.config = config or {}

    def is_ready(self):
        return False

    def summary(self):
        return {
            "name": "RNA-FM",
            "ready": False,
            "message": "Stub created. Implement model loader, preprocessing, inference, and validation."
        }
