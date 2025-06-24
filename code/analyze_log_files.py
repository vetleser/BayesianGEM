import logging
import re
import os
import csv
from collections import Counter, defaultdict
import pickle
from typing import Dict

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')
logging.info("BEGIN")

log_dir = "../results/analysis/apr14/"  # Change if needed
output_file = f"{log_dir}attempt_summary.csv"

logfile = "example_log.txt"

#
# Patterns

success_pattern = re.compile(r"Model solved successfully at temperature (\d{3}\.\d{2}) at attempt (\d+)")
failed_pattern = re.compile(r"Attempt 100 failed.*At Temperature (\d+\.\d{2})")


# Store attempts distribution: {temperature: Counter({attempt_count: frequency})}
attempts_distribution = defaultdict(list)
anaerobic_counter = 0

for number in range(100):
    logfile = f"{log_dir}gradient_search_{number}.log"
    condition = "Aerobic"
    try:
        with open(logfile, "r") as file:
            counter = 0
            for line in file:
                # Your custom logic goes here
                if counter >= 8:
                    condition = "Anaerobic"

                match1 = success_pattern.search(line)
                if match1:
                    counter += 1
                    temperature = float(match1.group(1))
                    attempts_count = int(match1.group(2))
                    attempts_distribution[(condition, temperature)].append(attempts_count)
                    if attempts_count > 10:
                        logging.info(f"Temperature {temperature:.2f} condition {condition} solved at attempt {attempts_count} in file {logfile}")
                    #logging.info(f"Temperature {temperature:.2f} solved at attempt {attempts_count}")
                match2 = failed_pattern.search(line)
                if match2:
                    counter += 1
                    temperature = float(match2.group(1))
                    attempts_count = 100
                    attempts_distribution[(condition, temperature)].append(attempts_count)

                    # logging.info(f"Temperature {temperature:.2f} condition {condition} failed at attempt {attempts_count} in file {logfile}")
                # For example, pattern matching or counting attempts
                

    except FileNotFoundError:
        logging.info(f"Log file not found: {logfile}")

# Log summary
logging.info("Temperature | Attempts to solve : Frequency")
for condition, attempts in attempts_distribution.items():
    number_of_attempts = len(attempts)
    if number_of_attempts != 100:
        logging.warning(f"Unexpected number of attempts for {condition[0]} at {condition[1]:.2f} K: {number_of_attempts}")
        continue
    attempts_counter = Counter(attempts)
    logging.info(f"{condition[0]} | {condition[1]:.2f} : {attempts_counter}")

