from .utils import Answer, log, LOG_ESSENTIAL


class Result:

    def __init__(self):
        self.PAST = Answer.UNKNOWN
        self.AST = Answer.UNKNOWN
        self.witnesses = []

    def all_known(self) -> bool:
        return self.PAST.is_known() and self.AST.is_known()

    def add_witness(self, witness):
        self.witnesses.append(witness)

    def print(self):
        log("", LOG_ESSENTIAL)
        log("", LOG_ESSENTIAL)
        log(f"PAST: {self.PAST}", LOG_ESSENTIAL)
        log(f"AST: {self.AST}", LOG_ESSENTIAL)
        log("", LOG_ESSENTIAL)
        log("", LOG_ESSENTIAL)
        for witness in self.witnesses:
            witness.print()
            log("", LOG_ESSENTIAL)
            log("", LOG_ESSENTIAL)