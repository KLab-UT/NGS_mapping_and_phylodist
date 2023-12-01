infile = open("mapped_percentage.txt", 'r')

class Sample:
    def __init__(self, name_in):
        self.name = name_in

    def add_total_merged(self, tmr_in):
        self.total_merged_reads = int(tmr_in)
       
    def add_total_unmerged(self, tur_in):
        self.total_unmerged_reads = int(tur_in)

    def add_mapped_merged(self, mmr_in):
        self.mapped_merged_reads = int(mmr_in)

    def add_mapped_unmerged(self, mur_in):
        self.mapped_unmerged_reads = int(mur_in)

    def get_combined_total(self):
        combined_total_reads = self.total_merged_reads + self.total_unmerged_reads
        return combined_total_reads
    
    def get_combined_mapped(self):
        combined_mapped_reads = self.mapped_merged_reads + self.mapped_unmerged_reads
        return combined_mapped_reads
    
    def add_merged_avgdepth(self, mad_in):
        self.merged_avgdepth = float(mad_in)

    def add_unmerged_avgdepth(self, uad_in):
        self.unmerged_avgdepth = float(uad_in)

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
        
        if ref_id not in ref_dict:
            ref_dict[ref_id] = {}
        if sample_id in ref_dict[ref_id]:
            # add to sample object in the dictionary
            if merge_type == "merged.bam":
                ref_dict[ref_id][sample_id].add_total_merged(total_reads)
                ref_dict[ref_id][sample_id].add_mapped_merged(mapped_reads)
                ref_dict[ref_id][sample_id].add_merged_avgdepth(avg_depth)
            elif merge_type == "unmerged.bam":
                ref_dict[ref_id][sample_id].add_total_unmerged(total_reads)
                ref_dict[ref_id][sample_id].add_mapped_unmerged(mapped_reads)
                ref_dict[ref_id][sample_id].add_unmerged_avgdepth(avg_depth)
        else:
            # create sample object and add it to the dictionary
            if merge_type == "merged.bam":
                ref_dict[ref_id][sample_id] = Sample(sample_id)
                ref_dict[ref_id][sample_id].add_total_merged(total_reads)
                ref_dict[ref_id][sample_id].add_mapped_merged(mapped_reads)
                ref_dict[ref_id][sample_id].add_merged_avgdepth(avg_depth)
            elif merge_type == "unmerged.bam":
                ref_dict[ref_id][sample_id] = Sample(sample_id)
                ref_dict[ref_id][sample_id].add_total_unmerged(total_reads)
                ref_dict[ref_id][sample_id].add_mapped_unmerged(mapped_reads)
                ref_dict[ref_id][sample_id].add_unmerged_avgdepth(avg_depth)


print(ref_dict)    
    
    
infile.close()
    
with open("tmp_outfile_RK.txt", 'w') as outfile:
    outfile.write(f"Reference_id,Sample_id,Total_reads,Mapped_reads,Percent_mapped,Avg_depth\n")
    for ref,sample_dict in ref_dict.items():
        for sample_id,sample_obj in sample_dict.items():
            percent_mapped = sample_obj.get_combined_mapped() / sample_obj.get_combined_total()
            outfile.write(f"{ref},{sample_id},{sample_obj.get_combined_total()},{sample_obj.get_combined_mapped()},{percent_mapped},{sample_obj.get_combined_avgdepth()}\n")
        
    
