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


IMAGE_CODE = 'asdasd3edwasd'


class CrypkoAttributes:
    def __init__(self, attr_string):
        self._attrs = int(attr_string)


class CrypkoOwner:
    def __init__(self, data):
        self.username = data.get('username')
        self.avatar = data.get('avatar')
        self.address = data.get('address')
        self.avatar_dat = data.get('avatarCrypko')
        self.bio = data.get('introduction')


class CrypkoDetails:
    def __init__(self, data, api):
        self.blush = data.pop('blush', False)
        self.dark_skin = data.pop('darkSkin', False)
        self.derivatives = data.pop('derivatives', [])
        self.hair_colour = data.pop('hairColor', '')[:-5]  # Strip " hair"
        self.hair_style = data.pop('hairStyle', '')
        self.hat = data.pop('hat', False)
        self.open_mouth = data.pop('openMouth', False)
        self.ribbon = data.pop('ribbon', False)
        self.smile = data.pop('smile', False)
        self.glasses = data.pop('glasses', False)

        self.updated_block_num = data.pop('updatedBlockNum', '')

        self.last_liked = data.pop('lastLiked', '')  # TODO: datetime

        self.sire = self.matron = None
        if 'sire' in data:
            self.sire = Crypko(data.pop('sire'), api)
        if 'matron' in data:
            self.matron = Crypko(data.pop('matron'), api)

        self.owner = CrypkoOwner(data.pop('owner'))

        self.bio = data.pop('bio', None)

        for i in ['id', 'attrs', 'cooldownIndex', 'eyeColor', 'iteration', 'likeCount',
                  'matronId', 'nextActionAt', 'noise', 'onRentalClockAuction',
                  'onSaleClockAuction', 'ownerAddr', 'cooldownRemain', 'sireId',
                  'name']:
            data.pop(i, None)

        # TODO: Handle these
        for i in ['deltaPriceDouble', 'duration', 'endAt', 'endingPrice', 'startedAt',
                  'startingPrice', 'startingPriceDouble']:
            data.pop(i, None)

        if data:
            print('[WARNING] Details has extra stuff: {}'.format(data))


class Crypko:
    def __init__(self, data, api):
        self._api = api

        self._details = None

        self.id = data.pop('id')
        self.attrs = CrypkoAttributes(data.pop('attrs', 0))

        self.likes = data.pop('likeCount', 0)
        self.noise = data.pop('noise', '')
        self.on_sale = data.pop('onSaleClockAuction', False)
        self.cooldown = data.get('cooldownRemain', 0) < 0
        self.duration = int(data.pop('duration', 0))
        self.end_price = int(data.pop('endingPrice', 0))
        self.iteration = data.pop('iteration', 0)
        self.on_rental = data.pop('onRentalClockAuction', False)
        self.eye_colour = data.pop('eyeColor', '')[:-5]  # Strip the " eyes"
        self.start_time = int(data.pop('startedAt', 0))
        self.start_price = int(data.pop('startingPrice', 0))
        self.cooldown_left = data.pop('cooldownRemain', 0)
        self.cooldown_index = data.pop('cooldownIndex')
        self.next_action_at = data.pop('nextActionAt', 0)

        self.name = data.pop('name', None)

        if data:
            print('[WARNING] Unhandled Crypko data: {}'.format(data))

    def __str__(self):
        rtn = 'Crypko #{}'.format(self.id)
        if self.name is not None:
            rtn += ' (' + self.name + ')'
        if self.on_sale:
            rtn += ' SALE'
        if self.on_rental:
            rtn += ' RENTAL'
        
        return rtn

    @property
    def details(self):
        if self._details is None:
            res = requests.get('{}crypkos/{}/detail'.format(self._api.DOMAIN, self.id))

            if res.status_code != 200:
                raise CrypkoApiError('Got response code {}. ({})'.format(res.status_code, res.text))

            self._details = CrypkoDetails(res.json(), self._api)

        return self._details

    @property
    def image(self):
        return 'https://img.crypko.ai/daisy/' + \
                hashlib.sha1((str(self.noise) + IMAGE_CODE + str(self.attrs._attrs)).encode()).hexdigest() + \
                '_lg.jpg'

    @property
    def image_sm(self):
        return self.image[:-6] + 'sm.jpg'

