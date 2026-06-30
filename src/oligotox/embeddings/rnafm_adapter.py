class RNAFMAdapter:
    def __init__(self):
        self.name = "rna_fm"

    def encode(self, sequences):
        raise NotImplementedError("Hook RNA-FM embedding extraction here.")
