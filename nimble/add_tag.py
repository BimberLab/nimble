import pysam
import sys

def add_concatenated_tag(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as infile:
        with pysam.AlignmentFile("-", "wb", template=infile) as outfile:
            for read in infile:
                try:
                    ur_tag = read.get_tag('UB')
                except KeyError:
                    try:
                        ur_tag = read.get_tag('UR')
                    except KeyError:
                        ur_tag = ''

                try:
                    cb_tag = read.get_tag('CB')
                except KeyError:
                    try:
                        cb_tag = read.get_tag('CR')
                    except KeyError:
                        cb_tag = ''

                concatenated_tag = ur_tag + cb_tag
                read.set_tag('XC', concatenated_tag)
                outfile.write(read)

if __name__ == "__main__":
    add_concatenated_tag(sys.argv[1])
