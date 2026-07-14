import pysam

from hcrseq.common.quantify import HCRseqQuantifier
from hcrseq.common.ref import Pathway, Reporter


class DummyReference:
    def __init__(self, reporters, pathways):
        self.reporters = reporters
        self.pathways = pathways


def build_dummy_reference():
    reporters = {
        name: Reporter(
            name=name,
            plasmid=name,
            sequence="A" * 10,
            tx_start=1,
            tx_end=10,
            barcode_position=1,
            barcode=None,
            lesion_position=1,
            unrepaired_base="A",
            repaired_base="G",
        )
        for name in ["control", "repairer"]
    }

    pathways = [
        Pathway(
            name="repair",
            metric="fraction_repaired",
            reporter=reporters["repairer"],
            control=reporters["control"],
        )
    ]

    return DummyReference(list(reporters.values()), pathways)


def test_quantify_repair_returns_pathway_measurements():
    reference = build_dummy_reference()
    quantifier = HCRseqQuantifier(reference, tags=None)

    quantifier.counters["control"].counts = {"total": 20, "repaired": 5, "unrepaired": 5}
    quantifier.counters["repairer"].counts = {"total": 10, "repaired": 6, "unrepaired": 4}

    measurements = quantifier.quantify_repair()

    assert set(measurements) == {"repair"}
    assert measurements["repair"]["value"] == 0.6
    assert measurements["repair"]["sd"] > 0


def test_count_umis_then_quantify_repair(tmp_path):
    reference = build_dummy_reference()
    quantifier = HCRseqQuantifier(reference, tags=None)

    bam_path = tmp_path / "synthetic.bam"
    fasta_path = tmp_path / "synthetic.fa"
    fasta_path.write_text(">control\nAAAAAA\n>repairer\nAAAAAA\n")

    header = {
        "HD": {"VN": "1.4"},
        "SQ": [{"SN": "control", "LN": 6}, {"SN": "repairer", "LN": 6}],
    }

    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
        for reporter_name in ["control", "repairer"]:
            for i in range(10):
                read = pysam.AlignedSegment()
                read.query_name = f"{reporter_name}_read_{i}"
                read.reference_id = bam_out.gettid(reporter_name)
                read.reference_start = 0
                read.cigar = [(0, 6)]
                if reporter_name == "repairer" and i < 6:
                    sequence = "G" + "A" * 5
                else:
                    sequence = "A" * 6
                read.query_sequence = sequence
                read.flag = 0
                read.mapq = 60
                read.tags = [("MD", f"{len(sequence)}")]
                bam_out.write(read)

    pysam.index(str(bam_path))

    quantifier.count_umis(str(bam_path))
    measurements = quantifier.quantify_repair()

    assert set(measurements) == {"repair"}
    assert measurements["repair"]["value"] == 0.6
    assert measurements["repair"]["sd"] > 0
