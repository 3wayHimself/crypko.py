import json
import time
import math

import crypko

from requests.adapters import HTTPAdapter


ADDR = '0xca39e90cec69838e73cc4f24ec5077dac44b47d6'
# ITER_WEIGHT = 1/10
# SIM_THRESHOLD = 12

CAP = 270000
USE_CAP = False

ATTRIBS = 'cooldownready'
FILTERS = None
# FILTERS = 'iteration:~3'  # Iter 3 or under only
# FILTERS = 'cooldown:ur'  # Only ultra rapids


def score(c1, c2):
    if USE_CAP and c1.id > CAP and c2.id > CAP:
        return 0

    if c1.on_sale or c2.on_sale:
        return 0

    similarity = bin(c1.attrs._attrs & c2.attrs._attrs).count('1')
    n_sim = bin(int(c1.noise) & int(c2.noise)).count('1')

    # if similarity < SIM_THRESHOLD:
    #     return 0
    # return similarity / (max(c1.iteration, c2.iteration) ** ITER_WEIGHT)

    return similarity * n_sim


def legal(c1, c2, cache):
    a = [c1.id, c2.id, *cache[c1.id], *cache[c2.id]]
    return len(a) == len(set(a))


def load_cache():
    try:
        with open('cache.json') as file_:
            dat = json.load(file_)
    except FileNotFoundError:
        dat = {}

    rtn = {}
    for i in dat:
        rtn[int(i)] = (*map(int, dat[i]),)

    return rtn


def save_cache(dat):
    with open('cache.json', 'w') as file_:
        json.dump(dat, file_)


def best_pair(crypkos, cooldown=False):
    scores = []
    cache = load_cache()

    for c1 in crypkos:
        if c1.id not in cache:
            cache[c1.id] = [c1.details.matron.id if c1.details.matron is not None else time.time(),
                    c1.details.sire.id if c1.details.sire is not None else time.time()]
            save_cache(cache)

        for c2 in crypkos:
            if c2.id not in cache:
                cache[c2.id] = [c2.details.matron.id if c2.details.matron is not None else time.time(),
                        c2.details.sire.id if c2.details.sire is not None else time.time()]
                save_cache(cache)

            if not legal(c1, c2, cache):
                continue

            scores.append((c1, c2, score(c1, c2)))
    
    scores.sort(key=lambda x:x[2], reverse=True)
    if scores:
        return scores[0][0], scores[0][1]
    return None


def purge_attributes(api, attribs):
    print(f'==> Collecting crypkos for \'{attribs}\'')
    crypkos = [i for i in api.search(owner_addr=ADDR, attributes=attribs)[1] if not i.on_sale]
    txs = []
    for crypko in crypkos:
        print(f'==> Selling #{crypko.id}')
        tx = api.contract.start_auc(crypko.id)
        print(f' ::  TX {tx.hex()}')
        txs.append(tx)

    print(f'==> Flushing {len(txs)} transactions')
    for tx in txs:
        print(f' ::  TX {tx.hex()}')
        api.contract.wait_for_tx(tx)
    print(f'==> Done!')


def main():
    with open('priv.key') as key_file:
        api = crypko.API(ADDR, key_file.read())

    batch_size = 12
    while True:
        purge_attributes(api, 'dark skin')
        purge_attributes(api, 'glasses')

        print(f'==> Requesting crypkos..')
        my_crypkos = [i for i in api.search(owner_addr=ADDR,
            attributes=ATTRIBS,
            filters=FILTERS,
            results=120)[1] if not i.on_sale]

        txs = []
        for _ in range(batch_size):
            print(f'==> Locating best pair..')
            pair = best_pair(my_crypkos)
            if pair is None:
                print(' ::  No legal fuses avaliable. Breaking.')
                break

            print(f' ::  Fusing #{pair[0].id} with #{pair[1].id}..')
            tx = api.contract.fuse(pair[0].id, pair[1].id)
            print(f' ::  Transaction {tx.hex()}')
            txs.append(tx)

            my_crypkos.remove(pair[0])
            my_crypkos.remove(pair[1])
        print(f'==> Flushing {len(txs)} transactions')
        for tx in txs:
            print(f' ::  TX {tx.hex()}')
            api.contract.wait_for_tx(tx)
        print(f' ::  Done!')


if __name__ == '__main__':
    main()

