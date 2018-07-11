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

import threading
import math

import requests

from .errors import CrypkoApiError, CrypkoNotFoundError
from .blockchain import ContractHandler
from .objects import Crypko


class API:
    DOMAIN = 'https://api.crypko.ai/'
    PER_PAGE = 12  # Can't change this :(

    def __init__(self, address=None, key=None, proxies=None):
        if key is not None and address is not None:
            self._contract = ContractHandler(key, address)
        else:
            self._contract = None
        self.proxies = proxies

    def _search(self, owner_addr=None, page=None, sort=None, category=None, attributes=None, filters=None):
        req = requests.get(self.DOMAIN + 'crypkos/search', params={
            'ownerAddr': owner_addr,
            'page': page,
            'sort': sort,
            'category': category,
            'attributes': attributes,
            'filters': filters
        }, proxies=self.proxies)
        if req.status_code != 200:
            raise CrypkoApiError('Got response code {}. ({})'.format(req.status_code, req.text))

        return req.json()

    @property
    def contract(self):
        if self._contract is None:
            raise CrypkoApiError('Contract not loaded. Did you provide an address and key?')
        return self._contract

    def get_crypko(self, crypko_id):
        res = requests.get('{}crypkos/{}/detail'.format(self.DOMAIN, crypko_id), proxies=self.proxies)

        if res.status_code == 404:
            raise CrypkoNotFoundError
        if res.status_code != 200:
            raise CrypkoApiError('Got response code {}. ({})'.format(res.status_code, res.text))

        return Crypko(res.json(), self)

    def threaded_search(self, callback, threads=8, start=0, results=-1, **kwargs):
        if results == 0:
            return
        start = max(0, start)
        page = math.floor(start / self.PER_PAGE) + 1

        kwargs.pop('page', None)
        if results < 0:
            r = self._search(**kwargs)
            results = r['totalMatched'] - start + 1

        page_end = math.ceil((start + results) / self.PER_PAGE) + 1

        def runner():
            nonlocal page
            while page < page_end:
                me_page = page
                page += 1
                res = self._search(page=me_page, **kwargs)

                if not res['crypkos']:
                    break
                
                for n, i in enumerate(res['crypkos']):
                    if n + (page - 1) * self.PER_PAGE >= start:
                        if n + (page - 1) * self.PER_PAGE < start + results:
                            callback(Crypko(i, self))

        pool = []
        for _ in range(threads):
            pool.append(threading.Thread(target=runner))
            pool[-1].start()

        for thread in pool:
            thread.join()

    def search(self, start=0, results=-1, **kwargs):
        """
        Search for Crypkos. Will retrieve `results` Crypkos starting at `start`.
        If results is negative, we will read until the end of results.

        TODO: Document kwargs
        """
        if results == 0:
            return
        start = max(0, start)
        page = math.floor(start / self.PER_PAGE) + 1
        kwargs.pop('page', None)

        result = self._search(page=page, **kwargs)

        def iters():
            nonlocal result, page
            results_sent = 0
            while results < 0 or results_sent < results:
                if not result['crypkos']:
                    break

                if results_sent == 0:
                    result['crypkos'] = result['crypkos'][start % self.PER_PAGE:]

                for i in result['crypkos']:
                    if results < 0 or results_sent < results:
                        yield Crypko(i, self)
                        results_sent += 1
                page += 1

                if results < 0 or results_sent < results:
                    result = self._search(page=page, **kwargs)
        
        return result['totalMatched'], iters()
