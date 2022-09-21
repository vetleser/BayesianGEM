#!/usr/bin/env python
# coding: utf-8
import logging
import pebble
import time
from pebble.concurrent.process import TimeoutError

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s') 

def my_fun(should_return):
    if should_return:
        return "Hei!"
    else:
        time.sleep(5)


logging.info("Initialize pool")
p = pebble.ProcessPool(2)
logging.info("Iterating")
raw_results = p.map(my_fun,[True, False,True,False,True],timeout=2).result()

results = []
while True:
    try:
        raw_res = raw_results.next()
    except StopIteration:
        break
    except TimeoutError:
        logging.info("A job timed out")
    else:
        logging.info("Job return successfully")
        results.append(raw_res)

print(results)
