import subprocess
import multiprocessing

def run_process(command, log_file):
    with open(log_file, "w") as f:
        print("hello")
        result = subprocess.run(command, shell=True, stdout=f, stderr=subprocess.STDOUT)

if __name__ == "__main__":
    # Define commands and log file names
    work_dir = "/data/rashika/CAFA4/"
    
    t0_gaf_file = work_dir + "uniprot/raw_goa/sample_t0.gz"
    #t0_gaf_file = work_dir + "uniprot/raw_goa/t0/goa_uniprot_all.gaf.195.gz"
    #t0_out_path = work_dir + "extracted_goa/t0_preprocessed.csv"
    t0_out_path = work_dir + "extracted_goa/t0_sample.csv"
    #log_t0 =  work_dir + "log/log_preprocess_t0.txt"
    log_t0 =  work_dir + "log/log_preprocess_t0_sample.txt"
    
    t1_gaf_file = work_dir + "uniprot/raw_goa/sample_t1.gz"
    #t1_gaf_file = work_dir + "uniprot/raw_goa/t1/goa_uniprot_all.gaf.gz"
    #t1_out_path = work_dir + "extracted_goa/t1_preprocessed.csv"
    t1_out_path = work_dir + "extracted_goa/t1_sample.csv"
    #log_t1 = work_dir + "log/log_preprocess_t1.txt"
    log_t1 =  work_dir + "log/log_preprocess_t1_sample.txt"
    
    cmd_preprocess_t0 = [
    "python3",                 # Command to execute Python 3
    "preprocess_gaf.py",       # Script to run
    t0_gaf_file,  # Path to input file
    "--highTP",
    "--out_path", t0_out_path,        # Output path parameter
    #"--evidence_codes", "EXP", "IDA",   # Evidence codes parameter
    #"--extract_col_list", "DB Object ID", "Qualifier"  # Extract column list parameter
    #" > ", log_t0,
]
    cmd_preprocess_t1 = [
    "python3",                 # Command to execute Python 3
    "preprocess_gaf.py",       # Script to run
    t1_gaf_file,  # Path to input file
    "--highTP",
    "--out_path", t1_out_path,        # Output path parameter
    #"--evidence_codes", "EXP", "IDA",   # Evidence codes parameter
    #"--extract_col_list", "DB Object ID", "Qualifier"  # Extract column list parameter
]
    

    # Create processes for each command
    process1 = multiprocessing.Process(target=run_process, args=(cmd_preprocess_t0, log_t0))
    process2 = multiprocessing.Process(target=run_process, args=(cmd_preprocess_t1, log_t1))
    #run_process(cmd_preprocess_t0, log_t0)
    #run_process(cmd_preprocess_t1, log_t0)
    # Start the processes
    process1.start()
    process2.start()

    # Wait for both processes to finish
    process1.join()
    process2.join()

    print("Both processes have finished.")
