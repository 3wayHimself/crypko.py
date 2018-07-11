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

# This example demonstrates how to use the threaded search to rapidly
# download the image for every crypko a user owns. It relies of a directory
# named out/ existing already.

import requests
import crypko

from requests.adapters import HTTPAdapter


THREADS = 16  # The number of seperate threads to run
# Replace this with the address of the user to download from
ADDRESS = '0xCa39E90CeC69838e73CC4F24Ec5077daC44B47d6'


def prepare_callback(session):
    def callback(cryp):
        print('Processing crypkoi#%d' % cryp.id)

        # cryp.image could be changes to cryp.image_sm to download the smaller
        # 192x192 version of the image instead of the 512x512 version.
        r = session.get(cryp.image, stream=True)

        # Download and save the image
        print(cryp.image)
        with open('out/%d.jpg' % cryp.id, 'wb') as file_:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    file_.write(chunk)

    return callback


def main():
    api = crypko.API()

    # Create a custom session so we can set a retry count
    session = requests.Session()
    session.mount('https://', HTTPAdapter(max_retries=5))

    # Perform the search
    api.threaded_search(prepare_callback(session), threads=THREADS, owner_addr=ADDRESS)


if __name__ == '__main__':
    main()

