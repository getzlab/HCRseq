import click
from hcrseq.amplicon.command_line import amplicon
from hcrseq.scRNA.command_line import scrna

@click.group()
def hcrseq():
    pass

hcrseq.add_command(amplicon)
hcrseq.add_command(scrna)