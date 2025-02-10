import sys
from collections import namedtuple



columns = ("seqname", "source", "feature",  "start", "end", "score", "strand", "frame", "attr",
            "rname", "rstart", "rend", "rscore")
out_cols = ("seqname", "source", "feature",  "start", "end", "score", "strand", "frame", "attr")
col_idx = dict(zip(columns, range(0, len(columns))))
class Record (namedtuple("Rec", columns)):
    @classmethod
    def from_line(cls, line):
        elts = line.strip().split('\t')
        for col in ("start", "end", "rstart", "rend"):
             elts[col_idx[col]] = int(elts[col_idx[col]])
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
        elif not rec.is_same_region(rec_prev):
                parse_region(rec_list, ofh)
                rec_prev = rec
                rec_list.clear()
                rec_list.append(rec)
        else:
            rec_list.append(rec)
    else:
        parse_region(rec_list, ofh)
    ifh.close()
    ofh.close()

def parse_region_alert(rec_list, ofh):
    score_list = list()
    rlength = 0
    # print(rec_list[0])
    try:
        if(int(rec_list[0].start) -1  < int(rec_list[0].rstart)):
             raise AssertionError("region start not match")
        if(rec_list[-1].end > rec_list[-1].end):
             raise AssertionError("region end not match")

        for i, rec in enumerate(rec_list):
            if i >= 1:
                if(rec.rstart != rec_list[i-1].rend):
                    raise AssertionError("blocks not connect to previous: %s" % rec.rstart)
            length = rec.rlength()
            score_list.extend([rec.rscore for i in range(max(rec.start -1, rec.rstart), min(rec.end, rec.rend))])

        if(rec_list[0].length() != len(score_list)):
             raise AssertionError("block length does not match")
        print(*(rec_list[0][col_idx[i]] for i in out_cols), score_list,
              sep = "\t", file=ofh)  
    # except OSError: 
    #    raise            
    except AssertionError as e:
         print(rec_list[0])
         print("AssersionError: %s" % str(e))

def parse_region(rec_list, ofh):
    known_score = dict()

    for i, rec in enumerate(rec_list):
        for pos in range(rec.rstart, rec.rend):
             known_score[pos] = rec.rscore
    score_list = [known_score.get(p, "-1") for p in range(rec_list[0].start -1 , rec_list[0].end)]
    score_str = ",".join(score_list)
    print(*(rec_list[0][col_idx[i]] for i in out_cols), score_str,
              sep = "\t", file=ofh)  



if __name__ == '__main__':
    if len(sys.argv) - 1 != 2:
        print("Usage: python <script> <in.merged.gtf> <out.with_score_list.gtf>")
        sys.exit(1)
    in_file, out_file = sys.argv[1:]
    parse_file(in_file, out_file)

    
        