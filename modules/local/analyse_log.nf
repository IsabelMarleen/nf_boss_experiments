// Analyse read log
process analyse_log {    
    input: 
        path log_file
    output: 
        path "*.csv", emit: analysed_log

    script:
    """
    analyse_boss_log.py $log_file ${log_file.getSimpleName()}_log_processed.csv
    """
}