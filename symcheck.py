# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 17:06:30 2019

@author: Roberto
"""

import argparse

from collections import defaultdict

refseq_names = 'name2'
hgnc_names = 'Approved symbol'


def parse_cmd():
    '''
    parsing command line parameters
    '''
    parser = argparse.ArgumentParser(description='A simple tool to verify gene names.')
    parser.add_argument('-v', '--version', action='version',
                        version='1.0.0')
    parser.add_argument('-r', '--refseq', type=str, required=True,
                        help='RefSeq data filename')
    parser.add_argument('-n', '--hgnc', type=str, required=True,
                        help='HGNC data filename')
    parser.add_argument('-f', '--file', type=str,
                        help='Text file containing gene symbols')
    parser.add_argument('-s', '--symbols', type=str,
                        help='String containing gene symbols')
    parser.add_argument('-m', '--mito', action='store_true',
                        help='Flag mitochondrial genes from HGNC')
    parser.add_argument('-a', '--ambiguous', action='store_true',
                        help='Show warnings for ambiguous gene symbols')
    parser.add_argument('-w', '--weird', action='store_true',
                        help='Flag alternative/unplaced chromosomes')
    parser.add_argument('-e', '--everything', action='store_true',
                        help='Show info for each gene name in the list')
    
    opts = parser.parse_args()
    return opts


# Helper functions

#TODO log the warnings instead of printing them
def log(*args):
    #print(*args)
    pass
        
def getinfo(linedata, info, columns, integer=False):
    item = linedata[columns[info]]
    if integer == False:
        return item
    else:
        return int(item)

def load_genes_from_file(filename):
    with open(filename, 'r', encoding='utf8') as f:
        data = [f.read().strip()]
    return separate(data, ['\n', '\t', ' '])
        
def load_genes_from_text(text):
    data = [text.strip()]
    return separate(data, ['+', '\t', ' '])
    
def separate(data, separators):
    for sep in separators:
        data = [y for x in data for y in x.split(sep)]
    data = [x.strip() for x in data]
    return data
   
    
def make_columns(header):
    return dict([(n, header.index(n)) for n in header])


def get_hgnc_aliases(linedata, hgnc_columns):
    alias = getinfo(linedata, 'Alias symbol', hgnc_columns)
    previous = getinfo(linedata, 'Previous symbol', hgnc_columns)
    return set([alias, previous])


class Gene:
    def __init__(self, linedata, normal_chromosomes, refseq_columns=None, hgnc_columns=None):
        self.hgnc = False
        self.refseq = False
        self.normal_chromosomes = normal_chromosomes
        # Establishing strings for __repr__
        self.txStart = '?'
        self.txEnd = '?'
        self.cyto = '?'
        
        if refseq_columns:
            self.load_refseq(linedata, refseq_columns)
        elif hgnc_columns:
            self.load_hgnc(linedata, hgnc_columns)
        else:
            raise ValueError("Either refseq_columns or hgnc_columns must be provided")
            


    def load_refseq(self, linedata, refseq_columns):
        if refseq_names not in refseq_columns:
            if hgnc_names in refseq_columns:
                raise ValueError("HGNC data provided instead of RefSeq data?")
            else:
                raise ValueError("Unknown format for refseq_columns.")
        self.refseq = True
        self.refseq_columns = refseq_columns
        self.symbol = getinfo(linedata, refseq_names, self.refseq_columns)
        self.chrom = getinfo(linedata, 'chrom', self.refseq_columns)
        self.refseq_update(linedata, refseq_columns, check=False)


    def load_hgnc(self, linedata, hgnc_columns):
        if hgnc_names not in hgnc_columns:
            if refseq_names in hgnc_columns:
                raise ValueError("RefSeq data provided instead of HGNC data?")
            else:
                raise ValueError("Unknown format for hgnc_columns.")
        self.hgnc = True
        self.hgnc_columns = hgnc_columns
        self.symbol = getinfo(linedata, hgnc_names, self.hgnc_columns)
        self.chrom = 'chr' + getinfo(linedata, 'Chromosome', self.hgnc_columns)
        self.hgnc_update(linedata, hgnc_columns, check=False)

        
    def refseq_update(self, linedata, refseq_columns, check=True):
        self.refseq_columns = refseq_columns
        dc = self.doublecheck(linedata, 'refseq')
        if check and not all(dc):
            if not dc[0] or self.chrom in self.normal_chromosomes \
                    or getinfo(linedata, 'chrom', self.refseq_columns) \
                    not in self.normal_chromosomes:
                log("WARNING: wrong merge attemped with gene", repr(self),
                    "and RefSeq line:\n"+str(linedata))
                return
            else: # This entry is placed in a conventional chromosome
                self.symbol = getinfo(linedata, 'chrom', self.refseq_columns)
                self.txStart = '?'
                self.txEnd = '?'
        self.refseq = True
        newstart = getinfo(linedata, 'txStart', self.refseq_columns, integer=True)
        newend = getinfo(linedata, 'txEnd', self.refseq_columns, integer=True)
        assert newstart < newend
        if self.txStart == '?' or self.txStart > newstart:
            self.txStart = newstart
        if self.txEnd == '?' or self.txEnd == None or self.txEnd < newend:
            self.txEnd = newend


    def hgnc_update(self, linedata, hgnc_columns, check=True):
        self.hgnc_columns = hgnc_columns
        dc = self.doublecheck(linedata, 'hgnc')
        if check and not all(dc):
            # I hate exceptions
            if not dc[0] or not (self.chrom in ['chrX', 'chrY'] and \
                                    getinfo(linedata, 'Chromosome', hgnc_columns) == 'X and Y'):
                log("WARNING: wrong merge attemped with gene", repr(self),
                    "and HGNC line:\n"+str(linedata))
                return
        self.hgnc = True
        self.status = getinfo(linedata, 'Status', self.hgnc_columns) #TODO remove this
        if self.cyto == '?':
            self.cyto = getinfo(linedata, 'Chromosome location', self.hgnc_columns)
        if not 'description' in self.__dir__():
            self.description = getinfo(linedata, 'Locus group', self.hgnc_columns)
        if not 'status' in self.__dir__():
            self.status = getinfo(linedata, 'Status', self.hgnc_columns)
        if not 'name' in self.__dir__():
            self.name = getinfo(linedata, 'Approved name', self.hgnc_columns)
        if not 'aliases' in self.__dir__():
            self.aliases = set()
        self.aliases.update(get_hgnc_aliases(linedata, self.hgnc_columns))



    def doublecheck(self, linedata, datatype):
        assert datatype in ['hgnc', 'refseq']
        if datatype == 'hgnc':
            symbol = getinfo(linedata, hgnc_names, self.hgnc_columns)
            chrom = 'chr' + getinfo(linedata, 'Chromosome', self.hgnc_columns)
        elif datatype == 'refseq':
            symbol = getinfo(linedata, refseq_names, self.refseq_columns)
            chrom = getinfo(linedata, 'chrom', self.refseq_columns)
        return [self.symbol == symbol, self.chrom == chrom]
                    
        
    def verbose(self):
        text = []
        text.append("Gene symbol: "+ self.symbol)
        text.append("Chromosome: " + self.chrom)
        text.append("RefSeq data: " + ["No", "Yes"][self.refseq])
        if self.refseq:
            text.append("Location: " + self.loc())
        text.append("HGNC data: " + ["No", "Yes"][self.hgnc])
        if self.hgnc:
            text.append("Cytoband: " + str(self.cyto))
            text.append("Status: " + str(self.status))
            text.append("Name: " + str(self.name))
            text.append("Locus group: " + str(self.description))
            text.append("Aliases of this specific gene: " + ', '.join(\
                        [item for item in self.aliases if item != '']))
        return text

    
    def loc(self):
        return '{}:{}-{}'.format(self.chrom, self.txStart, self.txEnd)


    def __repr__(self):
        return("<Gene:{}, {} {}>".format(self.symbol, self.loc(), self.cyto))



# Main functions

def make_databases(refseq_file, hgnc_file):
    print("Building database...")
    data = dict()
    refseq = set()
    hgnc_approved = set()
    hgnc_not_approved = set()
    aliases = defaultdict(set)
    normal_chromosomes = ['chr'+n for n in [str(i) for i in range(23)]+['X', 'Y']]
    for line in open(refseq_file, 'r', encoding="utf8"):
        if not line.startswith("#"):
            if "\t{}\t".format(refseq_names) in line:
                refseq_columns = make_columns(line.strip().split('\t'))
            else:
                linedata = line.strip().split('\t')
                symbol = getinfo(linedata, refseq_names, refseq_columns)
                if symbol in data:
                    data[symbol].refseq_update(linedata, refseq_columns)
                else:
                    refseq.add(symbol)
                    data[symbol] = Gene(linedata, normal_chromosomes,
                                        refseq_columns=refseq_columns)
    for i, line in enumerate(open(hgnc_file, 'r', encoding="utf8")):
        if not line.startswith("#"):
            if "\t{}\t".format(hgnc_names) in line:
                hgnc_columns = make_columns(line.strip().split('\t'))
            else:
                linedata = line.strip().split('\t')
                symbol = getinfo(linedata, hgnc_names, hgnc_columns)
                if getinfo(linedata, 'Status', hgnc_columns) != 'Approved':
                    hgnc_not_approved.add(symbol)
                else:
                    hgnc_approved.add(symbol)
                    for alias in [getinfo(linedata, 'Alias symbol', hgnc_columns),
                                  getinfo(linedata, 'Previous symbol', hgnc_columns)]:
                        if alias != '':
                            aliases[alias].add(symbol)
                    if symbol in data:
                        data[symbol].hgnc_update(linedata, hgnc_columns)
                    else:
                        data[symbol] = Gene(linedata, normal_chromosomes,
                                            hgnc_columns=hgnc_columns)
    return {'data': data,
            'hgnc_cols': hgnc_columns,
            'refseq_cols': refseq_columns, 
            'nl_chroms': normal_chromosomes,
            'refseq': refseq,
            'hgnc_approved': hgnc_approved,
            'hgnc_not_approved': hgnc_not_approved,
            'aliases': aliases}


def check_genelist(all_data):
    analysis = {"in_refseq": set(),
                "in_hgnc": set(), # can also be in refseq
                "withdrawn": set(),
                "refseq_missing": set(), # all missing from refseq
                "missing": set(), # missing in both refseq and HGNC
                "ambiguous": set(), # but not missing!
                "mito": set(),
                "weird": set()}
    for gene in set(all_data['genelist']):
        in_refseq = False
        missing = False
        if gene in all_data['refseq']:
            analysis['in_refseq'].add(gene)
            in_refseq = True
        else:
            analysis['refseq_missing'].add(gene)
        if gene in all_data['hgnc_approved']:
            analysis['in_hgnc'].add(gene)
        elif gene in all_data['hgnc_not_approved']:
            analysis['withdrawn'].add(gene)
        else:
            if not in_refseq:
                analysis['missing'].add(gene)
                missing = True
        if gene in all_data['aliases']:
            analysis['ambiguous'].add(gene)
        if not missing:
            if all_data['data'][gene].chrom in ['chrM', 'chrmitochondria']:
                analysis['mito'].add(gene)
            if all_data['data'][gene].chrom not in all_data['nl_chroms']:
                analysis['weird'].add(gene)
    return analysis
                
            
        
        

def check_repeated(genelist):
    return set([item for item in genelist if genelist.count(item) > 1])
        

    

def create_report(opts, all_data):
    print("Creating report...")
    with open("template", "r") as f:
        html = f.read()
    replaces = dict()
    replaces['refseq'] = st(all_data['refseq_file'])
    replaces['hgnc'] = st(all_data['hgnc_file'])
    replaces['genesno'] = len(all_data['genelist'])
    if opts.file:
        replaces['genesfile'] = "file: " + st(opts.file)
    if opts.symbols:
        replaces['genesfile'] = "string."

    rep = check_repeated(all_data['genelist'])
    if len(rep) >= 1:
        all_data['genelist'] = sorted(list(set(all_data['genelist'])))
        replaces['repwarn'] = "<br /><b>N.B.:</b> Please note that you have "\
                              "repeated the following names: " + ', '.join(rep) + \
                              "<br />The actual number of gene in your list is "\
                              "thus " + str(len(all_data['genelist'])) + "."
    else:
        replaces['repwarn'] = ''
    replaces['unique'] = len(all_data['genelist'])
    replaces['ambopt'] = [" do not", ""][bool(opts.ambiguous)]
    replaces['evopt'] = [" do not", ""][bool(opts.everything)]
    replaces['mitopt'] = [" <u>not</u>", ""][bool(opts.mito)]
    replaces['mitexpl'] = [" (you will not get any warning "\
                           "if a gene is mitochondrial)", ""][bool(opts.mito)]
    replaces['weirdopt'] = [" <u>not</u>", ""][bool(opts.weird)]
    replaces['weirdexpl'] = [" (you will not get any warning "\
                           "if a gene is in a non-standard chromosome)", ""][bool(opts.weird)]
    
    replaces['rsfound'] = len(all_data['analysis']['in_refseq'])
    replaces['hgncfound'] = len(all_data['analysis']['in_hgnc'].union(\
                                all_data['analysis']['withdrawn']).difference(\
                                all_data['analysis']['in_refseq']))
    if opts.mito:
        replaces['rsmito'] = " (" + str(len(all_data['analysis']['in_refseq'].intersection(\
                                all_data['analysis']['mito']))) + " as mitochondrial)"
        replaces['hgncmito'] = " (" + str(len((all_data['analysis']['in_hgnc'].union(\
                                all_data['analysis']['withdrawn'])).intersection(\
                                all_data['analysis']['mito']))) + " as mitochondrial)"
    else:
        replaces['rsmito'] = ''
        replaces['hgncmito'] = ''
    if opts.weird:
        replaces['rsweird'] = " (" + str(len(all_data['analysis']['in_refseq'].intersection(\
                                all_data['analysis']['weird']))) + " in non-standard chroms)"
        replaces['hgncweird'] = " (" + str(len((all_data['analysis']['in_hgnc'].union(\
                                all_data['analysis']['withdrawn'])).intersection(\
                                all_data['analysis']['weird']))) + " in non-standard chroms)"

    else:
        replaces['rsweird'] = ''
        replaces['hgncweird'] = ''

    replaces['rsmissingno'] = len(all_data['analysis']['refseq_missing'])
    replaces['hgncmissingno'] = len(all_data['analysis']['missing'])
    replaces['ambmissingno'] = len(all_data['analysis']['refseq_missing'].intersection(\
                                   all_data['analysis']['ambiguous']))
    
    if opts.ambiguous:
        replaces['ambtext'] = "Overall, " + str(len(all_data['analysis']['ambiguous'])) + \
                              " gene symbols were also aliases of other genes."
        amb = all_data['analysis']['ambiguous'].intersection(all_data['analysis']['in_refseq'])
        ambcore = 'a review of <b>all the symbols present in RefSeq but also having aliases</b>'
        replaces['ambiguous_review'] = possibly_empty_table(all_data,
                                                            amb,
                                                            ambcore)
    else:
        replaces['ambtext'] = ''
        replaces['ambiguous_review'] = ''

    if opts.mito:
        mit = all_data['analysis']['mito']
        mitcore = 'a list of <b>all mitochondrial genes (in "chrM" or "chrmitochondria" chromosomes)</b>'
        replaces['mitochondrial_review'] = possibly_empty_table(all_data,
                                                                mit,
                                                                mitcore)
    else:
        replaces['mitochondrial_review'] = ''
        
    if opts.weird:
        replaces['weirdtext'] = "Overall, " + str(len(all_data['analysis']['weird'])) + \
                              " genes were placed in non-standard chromosomes."
        wei = all_data['analysis']['weird']
        weicore = 'a review of <b>all the genes in non-standard chromosomes</b>'
        replaces['weird_review'] = possibly_empty_table(all_data,
                                                        wei,
                                                        weicore)
    else:
        replaces['weirdtext'] = ''
        replaces['weird_review'] = ''

    abs_ = list(set(all_data['genelist']).difference(all_data['analysis']['in_refseq']))
    replaces['absent_review'] = review(all_data, abs_)
    
    all_ = all_data['genelist']
    allt = '<p>Finally, by your request, here is a review of <b>all the unique symbols you provided:</b></p>'
    replaces['all_review'] = allt + review(all_data, all_)

    with open("output.html", "w") as f:
        f.write(html.format(**replaces))


def possibly_empty_table(all_data, geneset, coremessage):
    if len(geneset) > 0:
        text = '<p>As you asked, here is ' + coremessage + ':</b></p>'
        rev = review(all_data, geneset)
    else:
        text = '<p>You asked for ' + coremessage + '. There are none to display.'
        rev = ''
    return text + rev


def review(all_data, genelist):
    text = ['<table border="1">\n<tbody>\n']    
    genes = sorted(genelist)
    text.append("<tr>")
    text.append(td("No."))
    text.append(td("Your gene symbol"))
    text.append(td("RefSeq Location"))
    text.append(td("HGNC gene name"))
    text.append(td("HGNC locus group"))
    text.append(td("HGNC cytoband"))
    text.append(td("The symbol is also an alias of:"))

    for i, gene in enumerate(genes):
        text.append("<tr>")
        text.append(td(i+1))
        text.append(td(tooltip(all_data, gene, where='right')))
        if gene in all_data['analysis']['in_refseq']:
            text.append(td(all_data['data'][gene].loc()))
        else:
            text.append(td("absent"))
        if gene in all_data['analysis']['in_hgnc']:
            text.append(td(all_data['data'][gene].name))
            text.append(td(all_data['data'][gene].description))
            text.append(td(all_data['data'][gene].cyto))
        elif gene in all_data['analysis']['withdrawn']:
            text.append(td("(Withdrawn)"))
            text.append(td("-"))
            text.append(td("-"))
        else:
            text.append(td("(Missing)"))
            text.append(td("-"))
            text.append(td("-"))
        if gene in all_data['analysis']['ambiguous']:
            ambtool = [tooltip(all_data, gene, where='left') for gene in all_data['aliases'][gene]]
            text.append(td(', '.join(ambtool)))
        else:
            text.append(td("-"))
        text.append("</tr>")
    text.append("</tbody>\n</table>\n")
    return '\n'.join(text)

   
def td(text):
    return "<td><div>"+str(text)+"</div></td>"


def add_ghr_link(genename):
    return '<a href="https://ghr.nlm.nih.gov/search?query='+genename+'" target="_blank">'+genename+'</a>'


def tooltip(all_data, genename, where='left', add_link=True):
    if where == 'right':
        text = ['<div class="rtooltip">']
    else:
        text = ['<div class="ltooltip">']
    link = add_ghr_link(genename)
    text.append(link)
    text.append('<span class="tooltiptext"><p>')
    if genename in all_data['data'].keys():
        for line in all_data['data'][genename].verbose():
            text.append("{}<br />".format(line))
    else:
        text.append("(Not found)")
    text.append('</p></span></div>')
    return (''.join(text)) 


def st(text):
    '''
    sanitizes text
    '''
    return text.replace("<", "\<")

def main():
    opts = parse_cmd()
    refseq_file = opts.refseq
    hgnc_file = opts.hgnc
    if opts.file:
        genelist = load_genes_from_file(opts.file)
    elif opts.symbols:
        genelist = load_genes_from_text(opts.symbols)
    else:
        print("You need to specify gene names either from a file (-f) or as string (-s)")
        return
    print("\nUsing refseq database:", refseq_file.split('/')[-1])
    print("Using hgnc database:", hgnc_file.split('/')[-1])
    all_data = make_databases(refseq_file, hgnc_file)
    all_data['refseq_file'] = refseq_file
    all_data['hgnc_file'] = hgnc_file
    all_data['genelist'] = sorted(genelist)
    all_data['analysis'] = check_genelist(all_data)
    create_report(opts, all_data)
    return all_data
    


if __name__=='__main__':
    all_data = main()
