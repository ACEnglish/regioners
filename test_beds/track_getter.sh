#Downloads to beside this script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# EPD promoters
python -c """
import urllib.request
import json
import pandas as pd

url = 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=epdNewPromoter'
with urllib.request.urlopen(url) as response:
    html = response.read()
    data = json.loads(html)['epdNewPromoter']
    for i in data:
        print(i['chrom'], i['chromStart'], i['chromEnd'], sep='\t')
""" > $DIR/grch38.epd_promoters.bed

# CpG islands
python -c """
import urllib.request
import json
import itertools

url = 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cpgIslandExt'
with urllib.request.urlopen(url) as response:
    html = response.read()
    data = itertools.chain(*json.loads(html.decode())['cpgIslandExt'].values())
    for i in data:
        if '_' not in i['chrom']:
            print(i['chrom'], i['chromStart'], i['chromEnd'], sep='\t')
""" > $DIR/grch38.cpg_islands.bed
