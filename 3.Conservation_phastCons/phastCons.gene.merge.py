import sys
from collections import namedtuple



columns = ("seqname", "source", "feature",  "start", "end", "score", "strand", "frame", "attr",
            "rscore")
#out_cols = ("seqname", "source", "feature",  "start", "end", "score", "strand", "frame", "attr")
col_idx = dict(zip(columns, range(0, len(columns))))
class Record (namedtuple("Rec", columns)):
    @classmethod
    def from_line(cls, line):
        elts = line.strip().split('\t')
        for col in ("start", "end"):
             elts[col_idx[col]] = int(elts[col_idx[col]])
        elts[col_idx['rscore']] = str.split(elts[col_idx['rscore']], ",")
        return cls(*elts)
         
    def length(self):
        return self.end - self.start + 1
    def rlength(self):
        return self.rend - self.rstart
    def is_same_region(self, other):
        return all((self.attr == other.attr,
                    self.start == other.start,
                    self.end == other.end))

def parse_file(in_file, out_file):
    # file = "./ensembl_canonical.cds.phastCons.gtf"
    ifh = open(in_file, 'r')
    ofh = open(out_file, 'w')
    rec_prev = None
    rec_list = list()
    for line in ifh:
        rec = Record.from_line(line)
        if rec_prev is None:
                rec_prev = rec
                rec_list.append(rec)
                continue
        elif not rec.attr == rec_prev.attr:
                parse_gene(rec_list, ofh)
                rec_prev = rec
                rec_list.clear()
                rec_list.append(rec)
        else:
            rec_list.append(rec)
    else:
        parse_gene(rec_list, ofh)
    ifh.close()
    ofh.close()



def parse_gene(rec_list, ofh):
    #known_score = dict()
    rscore = list()
    length = 0
    strand = rec_list[0].strand
    gene = rec_list[0].attr

    for i, rec in enumerate(rec_list):
        length += rec.length()
        rscore.extend(rec.rscore)
        if (rec.strand != strand):
             print("%s: have both strand" % gene)
    if (length % 3 != 0):
         print("%s : length is not multiple of 3: %d" % (gene, length))
    if (length != len(rscore)):
         print("%s: gene length and score length not match" % gene)
    if strand == '-':
         rscore = rscore[::-1]
    rscore_str = ",".join(rscore)
    print(gene, length, rscore_str,
              sep = "\t", file=ofh)  



if __name__ == '__main__':
    if len(sys.argv) - 1 != 2:
        print("Usage: python <script> <cds.score_list.sort_by_gene.gtf> <out.gene.score_list.gtf>")
        sys.exit(1)
    in_file, out_file = sys.argv[1:]
    parse_file(in_file, out_file)

    
        