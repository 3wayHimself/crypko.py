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

import requests
import hashlib

from .errors import CrypkoApiError

IMAGE_CODE = 'asdasd3edwasd'

ATTR_MAP = {'hair': ['blonde hair', 'brown hair', 'black hair', 'blue hair', 'pink hair', 'purple hair', 'green hair',
                     'red hair', 'silver hair', 'white hair', 'orange hair', 'aqua hair', 'grey hair'],
            'eyes': ['blue eyes', 'red eyes', 'brown eyes', 'green eyes', 'purple eyes', 'yellow eyes', 'pink eyes',
                     'aqua eyes', 'black eyes', 'orange eyes'],
            'hair_style': ['long hair', 'short hair', 'twintails', 'drill hair', 'ponytail'],
            'attribute': ['dark skin', 'blush', 'smile', 'open mouth', 'hat', 'ribbon', 'glasses']}


class CrypkoAttributes:
    def __init__(self, attr_string):
        self.attribute_int = int(attr_string)

    @property
    def attributes(self):
        attrs = bin(self.attribute_int)[2:].rjust(38, '0')

        rtn = [
            (ATTR_MAP['hair'][int(attrs[30: 34], 2) % len(ATTR_MAP['hair'])], True),
            (ATTR_MAP['eyes'][int(attrs[22: 26], 2) % len(ATTR_MAP['eyes'])], True),
            (ATTR_MAP['hair_style'][int(attrs[14: 18], 2) % len(ATTR_MAP['hair_style'])], True)
        ]

        for n, i in enumerate(ATTR_MAP['attribute'][::-1]):
            rtn.append(
                (i, bool(int(attrs[n * 2])))
            )

        return rtn

    @property
    def sub_attributes(self):
        attrs = bin(self.attribute_int)[2:].rjust(38, '0')

        rtn = [
            (ATTR_MAP['hair'][int(attrs[34: 38], 2) % len(ATTR_MAP['hair'])], True),
            (ATTR_MAP['eyes'][int(attrs[26: 30], 2) % len(ATTR_MAP['eyes'])], True),
            (ATTR_MAP['hair_style'][int(attrs[18: 22], 2) % len(ATTR_MAP['hair_style'])], True)
        ]

        for n, i in enumerate(ATTR_MAP['attribute'][::-1]):
            rtn.append(
                (i, bool(int(attrs[n * 2 + 1])))
            )

        return rtn


class CrypkoOwner:
    def __init__(self, data):
        self.avatar_dat = data.get('avatarCrypko')
        self.username = data.get('username')
        self.address = data.get('address')
        self.avatar = data.get('avatar')
        self.bio = data.get('introduction')


class CrypkoAuction:
    SALE = 0
    RENTAL = 1

    def __init__(self, start_price, end_price, duration, started, ends, sale, rental):
        self.start_price = start_price
        self.end_price = end_price
        self.duration = duration
        self.started = started
        self.ends = ends

        self.active = sale or rental

        self.type = self.SALE if sale else self.RENTAL if rental else None


class Crypko:
    def __init__(self, data, api):
        self._api = api

        self.sire = self.matron = self.owner = None

        self.complete = False

        self._load_details(data)

    def _load_details(self, data):
        self.id = data.pop('id')

        # Attribute parsing
        attrs = CrypkoAttributes(data.pop('attrs', 0))
        self.sub_attributes = attrs.sub_attributes
        self.attributes_int = attrs.attribute_int
        self.attributes = attrs.attributes

        # Load auction details
        start_price = int(data.pop('startingPrice', 0))
        end_price = int(data.pop('endingPrice', 0))
        duration = int(data.pop('duration', 0))
        started = int(data.pop('startedAt', 0))
        ends = int(data.pop('endAt', 0))

        on_sale = data.pop('onSaleClockAuction', False)
        on_rental = data.pop('onRentalClockAuction', False)

        self.auction = CrypkoAuction(start_price, end_price, duration, started, ends, on_sale, on_rental)

        # Metadata
        self.derivatives = data.pop('derivatives', [])
        self.likes = data.pop('likeCount', 0)
        self.raw_noise = int(data.pop('noise', 0))

        self.iteration = data.pop('iteration', 0)

        self.cooldown_left = data.pop('cooldownRemain', 0)
        self.cooldown_index = data.pop('cooldownIndex')

        if 'owner' in data:
            self.owner = CrypkoOwner(data.pop('owner'))
            data.pop('ownerAddr', None)
        elif 'ownerAddr' in data:
            self.owner = CrypkoOwner({'address': data.pop('ownerAddr')})

        self.bio = data.pop('bio', None)
        self.name = data.pop('name', None)

        self.next_action_at = data.pop('nextActionAt', 0)

        # Visual attributes
        # # Features
        self.hair_colour = data.pop('hairColor', '')[:-5] or None
        self.eye_colour = data.pop('eyeColor', '')[:-5] or None
        self.hair_style = data.pop('hairStyle', None)
        self.open_mouth = data.pop('openMouth', None)
        self.dark_skin = data.pop('darkSkin', None)
        self.blush = data.pop('blush', None)
        # # Objects
        self.glasses = data.pop('glasses', None)
        self.ribbon = data.pop('ribbon', None)
        self.smile = data.pop('smile', None)
        self.hat = data.pop('hat', None)

        self.updated_block_num = data.pop('updatedBlockNum', '')

        self.last_liked = data.pop('lastLiked', '')  # TODO: datetime

        # Genetics
        if 'sire' in data:
            self.sire = Crypko(data.pop('sire'), self._api)
            data.pop('sireId', None)
        elif 'sireId' in data:
            self.sire = Crypko({'id': data.pop('sireId')}, self._api)

        if 'matron' in data:
            self.matron = Crypko(data.pop('matron'), self._api)
            data.pop('matronId', None)
        elif 'matronId' in data:
            self.matron = Crypko({'id': data.pop('matronId')}, self._api)

        if data:
            print('[WARNING] Unhandled Crypko data: {}'.format(data))

    def __str__(self):
        rtn = 'Crypko #{}'.format(self.id)
        if self.name is not None:
            rtn += ' (' + self.name + ')'
        if self.auction.active:
            if self.auction.type == CrypkoAuction.SALE:
                rtn += ' (on sale)'
            else:
                rtn += ' (for rent)'

        return rtn

    def __eq__(self, other):
        return other.id == self.id

    @property
    def noise(self):
        noise = bin(self.raw_noise)[2:].rjust(256, '0')
        return [int(noise[n - 2: n], 2) for n in range(len(noise), 0, -2)]

    def reload_details(self):
        res = requests.get('{}crypkos/{}/detail'.format(self._api.DOMAIN, self.id))

        if res.status_code != 200:
            raise CrypkoApiError('Got response code {}. ({})'.format(res.status_code, res.text))

        self._load_details(res.json())
        self.complete = True

    def ensure_complete(self):
        if not self.complete:
            self.reload_details()

    @property
    def image(self):
        return 'https://img.crypko.ai/daisy/' + \
               hashlib.sha1((str(self.noise) + IMAGE_CODE + str(self.attrs.attribute_int)).encode()).hexdigest() + \
               '_lg.jpg'

    @property
    def image_sm(self):
        return self.image[:-6] + 'sm.jpg'
