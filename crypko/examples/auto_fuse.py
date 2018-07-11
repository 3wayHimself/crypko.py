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

score_cache = {}

# FILTERS = 'iteration:~3'  # Iter 3 or under only
# FILTERS = 'cooldown:ur'  # Only ultra rapids


def score(c1, c2):
    if USE_CAP and c1.id > CAP and c2.id > CAP:
        return 0

    if c1.auction.active or c2.auction.active:
        return 0

    if (c1.id, c2.id) in score_cache:
        return score_cache[(c1.id, c2.id)]

    similarity = len((*filter(c1.attributes.__contains__, c2.attributes),))

    n_sim = sum(x == y for x, y in zip(c1.noise, c2.noise))

    # if similarity < SIM_THRESHOLD:
    #     return 0
    # return similarity / (max(c1.iteration, c2.iteration) ** ITER_WEIGHT)

    score_cache[(c1.id, c2.id)] = similarity * n_sim
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
            c1.ensure_complete()
            cache[c1.id] = (c1.matron.id if c1.matron is not None else time.time(),
                            c1.sire.id if c1.sire is not None else time.time())
            save_cache(cache)

        for c2 in crypkos:
            if c2.id not in cache:
                c2.ensure_complete()
                cache[c2.id] = (c2.matron.id if c2.matron is not None else time.time(),
                                c2.sire.id if c2.sire is not None else time.time())
                save_cache(cache)

            if not legal(c1, c2, cache):
                continue

            scores.append((c1, c2, score(c1, c2)))

    scores.sort(key=lambda x: x[2], reverse=True)

    if scores:
        return scores[0][0], scores[0][1]
    return None


def purge_attributes(api, attribs):
    print(f'==> Collecting crypkos for \'{attribs}\'')
    crypkos = [i for i in api.search(owner_addr=ADDR, attributes=attribs)[1] if not i.auction.active]
    txs = []
    for c in crypkos:
        print(f'==> Selling #{c.id}')
        tx = c.sell(0.7, 0.4)
        print(f' ::  TX {tx.hex()}')
        txs.append(tx)

    print(f'==> Flushing {len(txs)} transactions')
    for tx in txs:
        print(f' ::  TX {tx.hex()}')
        api.contract.wait_for_tx(tx)
    print(f'==> Done!')


def main():
    with open('priv.key') as key_file:
        api = crypko.API(ADDR, key_file.read().strip())

    batch_size = 48
    while True:
        purge_attributes(api, 'dark skin')
        purge_attributes(api, 'glasses')

        print(f'==> Requesting crypkos..')
        my_crypkos = [i for i in api.search(owner_addr=ADDR,
                                            attributes=ATTRIBS,
                                            results=12 * 300,
                                            filters=FILTERS)[1] if not i.auction.active]

        print(f' ::  {len(my_crypkos)} usable')

        txs = []
        for _ in range(batch_size):
            print(f'==> Locating best pair..')
            pair = best_pair(my_crypkos)
            if pair is None:
                print(' ::  No legal fuses avaliable. Breaking.')
                break

            print(f' ::  Fusing #{pair[0].id} with #{pair[1].id}..')
            try:
                tx = pair[0].fuse(pair[1])
                print(f' ::  Transaction {tx.hex()}')
                txs.append(tx)

                my_crypkos.remove(pair[0])
                my_crypkos.remove(pair[1])
            except ValueError as e:
                print(f' ::  Failed with {e}')
        print(f'==> Flushing {len(txs)} transactions')
        for tx in txs:
            print(f' ::  TX {tx.hex()}')
            api.contract.wait_for_tx(tx)
        print(f' ::  Done!')


if __name__ == '__main__':
    main()
