infile = open("mapped_percentage.txt", 'r')

class Sample:
    __init__(self, name_in):
        self.name = name_in

    add_total_merged(self, tmr_in):
        self.total_merged_reads = tmr_in
       
    add_total_unmerged(self, tumr_in):
        self.total_unmerged_reads = tumr_in

# add methods to add mapped_merged and mapped_unmerged, and then to get combined proportion


ref_dict = {}

for line in infile:
    if line[0] != "#":
        line_split = line.split(",")
        sample_id = line_split[0]
        ref_id = line_split[1]
        merge_type = line_split[2]
        total_reads = line_split[3]
        mapped_reads = line_split[4]
        if ref_id in ref_dict:
            if sample_id in ref_dict[ref_id]:
                # adding to the sample object that already exists
            else:
                # create sample object and add it to the dictionary


infile.close()
