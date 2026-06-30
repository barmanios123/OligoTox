class RiNALMoAdapter:
    def __init__(self):
        self.name = "rinalmo"

    def encode(self, sequences):
        raise NotImplementedError("Hook RiNALMo embedding extraction here.")
