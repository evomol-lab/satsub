from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def make_alignment(seqs: dict[str, str]) -> MultipleSeqAlignment:
    records = [SeqRecord(Seq(s), id=name, name=name, description="") for name, s in seqs.items()]
    return MultipleSeqAlignment(records)
