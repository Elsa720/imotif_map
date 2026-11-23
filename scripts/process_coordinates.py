import os

def process_narrowpeak(input_file, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            name = parts[3] if len(parts) > 3 else '.'
            score = parts[4] if len(parts) > 4 else '.'
            strand = parts[5] if len(parts) > 5 else '.'
            summit_offset = int(parts[9])
            
            summit = start + summit_offset
            new_start = max(0, summit - 250)
            new_end = summit + 250
            
            outfile.write(f"{chrom}\t{new_start}\t{new_end}\t{name}\t{score}\t{strand}\n")

def process_bed(input_file, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            center = (start + end) // 2
            new_start = max(0, center - 250)
            new_end = center + 250
            
            outfile.write(f"{chrom}\t{new_start}\t{new_end}\n")

if __name__ == "__main__":
    peak_in = "/data1/minj/footprint/GM_stast/basic_feature/0.rawdata/iMab.idr0.05.narrowPeak"
    peak_out = "/data1/minj/footprint/GM_stast/basic_feature/0.rawdata/meme/iMab.idr0.05.summit_500bp.bed"
    
    bg_in = "/data1/minj/footprint/GM_stast/basic_feature/0.rawdata/bg/bg.unified_matched.clean.bed"
    bg_out = "/data1/minj/footprint/GM_stast/basic_feature/0.rawdata/meme/bg.unified_matched.clean.center_500bp.bed"
    
    print("Processing narrowPeak file...")
    process_narrowpeak(peak_in, peak_out)
    
    print("Processing background bed file...")
    process_bed(bg_in, bg_out)
    
    print("Files successfully written to /data1/minj/footprint/GM_stast/basic_feature/0.rawdata/meme")
