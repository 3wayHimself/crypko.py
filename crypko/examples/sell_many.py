import crypko


ADDR = '0xca39e90cec69838e73cc4f24ec5077dac44b47d6'
SELL_DOWN_TO = 251824

GIFT = True

CAN_FIRST = 100
GIVE_TO = '0xb24978d166e766d9557fa75d9b7b16f0beccebc2'

START_PRICE = 0.8
END_PRICE = 0.6


def main():
    with open('priv.key') as key_file:
        api = crypko.API(ADDR, key_file.read().strip())

    txs = []
    count, crypkos = api.search(owner_addr=ADDR)

    cnt = 0
    for c in crypkos:
        if GIFT and not c.auction.active:
            if c.id not in [366013] and c.id > SELL_DOWN_TO:
                print(f'==> Gifting #{c.id}')

                tx = api.contract.gift(c.id, GIVE_TO)
                print(f' ::  TX {tx.hex()}')
                txs.append(tx)

                cnt += 1

                if cnt >= CAN_FIRST:
                    break
        elif not GIFT and c.id > SELL_DOWN_TO and not c.on_sale:
            print(f'==> Selling #{c.id}')

            tx = api.contract.start_auc(c.id, START_PRICE, END_PRICE)
            print(f' ::  TX {tx.hex()}')
            txs.append(tx)

    print(f'==> Flushing {len(txs)} transactions')
    for tx in txs:
        print(f' ::  TX {tx.hex()}')
        api.contract.wait_for_tx(tx)

    print(f'==> Done!')


if __name__ == '__main__':
    main()
