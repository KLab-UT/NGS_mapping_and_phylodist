infile = open("mapped_percentage.txt", 'r')

class Sample:
    def __init__(self, name_in):
        self.name = name_in

    def add_total_merged(self, tmr_in):
        self.total_merged_reads = tmr_in
       
    def add_total_unmerged(self, tur_in):
        self.total_unmerged_reads = tur_in

    def add_mapped_merged(self, mmr_in):
        self.mapped_merged_reads = mmr_in

    def add_mapped_unmerged(self, mur_in):
        self.mapped_unmerged_reads = mur_in

    def get_combined_total(self):
        combined_total_reads = self.total_merged_reads + self.total_unmerged_reads
        return combined_total_reads
    
    def get_combined_mapped(self):
        combined_mapped_reads = self.mapped_merged_reads + self.mapped_unmerged_reads
        return combined_mapped_reads
    
    def add_merged_avgdepth(self, mad_in):
        self.merged_avgdepth = mad_in

    def add_unmerged_avgdepth(self, uad_in):
        self.unmerged_avgdepth = uad_in

    def get_combined_avgdepth(self):
        combined_avgdepth = self.merged_avgdepth + self.unmerged_avgdepth
        return combined_avgdepth

# add methods to add mapped_merged and mapped_unmerged, and then to get combined proportion


ref_dict = {}

for line in infile:
    if line[0] != "#":
        line_split = line.split(",")
        sample_id = line_split[0]
        ref_id = line_split[1]
        merge_type = line_split[2]
        total_reads = line_split[3]
        avg_depth = line_split[4]
        mapped_reads = line_split[5]
        sample = Sample(sample_id)
        if merge_type = "merged.bam":
            sample.add_total_merged(total_reads)
            sample.add_mapped_merged(mapped_reads)
            sample.add_merged_avgdepth(avg_depth)
        elif merge_type = "unmerged.bam":
            sample.add_total_merged(total_reads)
            sample.add_mapped_merged(mapped_reads)
            sample.add_merged_avgdepth(avg_depth)
        
        if ref_id not in ref_dict:
            ref_dict[ref_id] = {}
        if sample_id in ref_dict[ref_id]:
            # adding to the sample object that already exists
            pass
        else:
            # create sample object and add it to the dictionary
            ref_dict[ref_id][sample_id] = Sample(sample_id)

print(ref_dict)    
    
    
infile.close()
    
    
