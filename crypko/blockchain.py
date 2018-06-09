# Copyright (c) 2018 Bottersnike
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import json
import os

from web3.middleware import geth_poa_middleware
from web3 import Web3, HTTPProvider


class ContractHandler:
    NODE = 'https://mainnet.infura.io/metamask'
    NODE = 'https://rinkeby.infura.io/metamask'  # Remove when switched to mainnet
    CONTRACT_ADDRESS = '0xFe2149773B3513703E79Ad23D05A778A185016ee'

    SPECIAL_ADDRS = {
        'ClockAuction': '0xC8cD5cA4f1a91cdb103EC19d3B46a1E42A862025'
    }

    RECEPIENT = '0xfe2149773b3513703e79ad23d05a778a185016ee'

    def __init__(self, key, address):
        self.web3 = Web3(HTTPProvider(self.NODE))
        self.web3.middleware_stack.inject(geth_poa_middleware, layer=0)

        # Load the ABI
        abi_path = os.path.join(os.path.dirname(__file__), 'abi.json')
        with open(abi_path) as file_:
            abi = json.load(file_)

        self.contracts = {}
        for i in abi:
            self.contracts[i] = self.web3.eth.contract(
                    address=self.SPECIAL_ADDRS.get(i, self.CONTRACT_ADDRESS),
                    abi=abi[i])

        # Variables for later
        self.priv_key = key
        self.address = address
        self.nonce = None

    def _send(self, abi, value=0):
        if self.nonce is None:
            self.nonce = self.web3.eth.getTransactionCount(Web3.toChecksumAddress(self.address))
        else:
            self.nonce += 1

        tx = { 
            'nonce': self.nonce,
            'gasPrice': Web3.toWei(1, 'gwei'),
            'from': Web3.toChecksumAddress(self.address),
            'gas': 1000000000,
            'value': value,
        }   
        gas = abi.estimateGas(tx)
        tx['gas'] = int(gas * 1.3)
        tx = abi.buildTransaction(tx)
        tx = self.web3.eth.account.signTransaction(tx, self.priv_key)

        return self.web3.eth.sendRawTransaction(tx.rawTransaction)

    def start_auc(self, crypko, start=0.7, end=0.4, duration=172800):
        start = Web3.toWei(start, 'ether')
        end = Web3.toWei(end, 'ether')
        start, end = int(start), int(end)

        create = self.contracts['CrypkoCore'].functions.createSaleAuction(crypko, start, end, duration)

        return self._send(create)

    def bid(self, crypko):
        cost = self.contracts['ClockAuction'].getCurrentPrice(crypko).call()
        bid = self.contracts['ClockAuction'].functions.bid(crypko)

        return self._send(bid, cost)

    def fuse(self, first, second):
        fee = self.contracts['CrypkoFusion'].call().fusionFee()
        fwa = self.contracts['CrypkoFusion'].functions.fuseWithAuto(first, second)

        return self._send(fwa, fee)

    def wait_for_tx(self, tx):
        return self.web3.eth.waitForTransactionReceipt(tx)
