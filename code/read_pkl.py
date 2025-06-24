import pickle
import logging
import os
from reframed import CBModel

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

path = os.path.dirname(os.path.realpath(__file__)).replace('code','')

model: CBModel = pickle.load(open(os.path.join(path, 'models/anaerobic.pkl'), 'rb'))

for rxn in model.reactions.values():
    if "oxygen" in rxn.name.lower():
        CBModel.print_reaction(model, rxn.id, use_names=True)
