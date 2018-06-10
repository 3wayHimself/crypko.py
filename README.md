# crypko.py

Install using `pip install crypko.py`.

Provides a simple wrapper around the Crypko platform allowing users to
search for crypko, purchase crypko, fuse crypko, etc.

More features are slowly being added as I work on them. If any part of
this stops working because the Crypko platform has been updated, let
me know and I'll try and fix it.

## Basic usage:
```py
>>> import crypko
>>> api = crypko.API()
>>> num, results = api.search(owner_addr='0xCa39E90CeC69838e73CC4F24Ec5077daC44B47d6', attributes='glasses')
>>> num
116
>>> next(resuls)
<crypko.objects.Crypko object at 0x7f01d9ded240>
>>> print(next(results))
Crypko #168993 SALE
>>> c = next(results)
>>> c.image
'https://img.crypko.ai/daisy/9753747176b22ea4fdb120f308073c35930c38a6_lg.jpg'
>>> c.details.matron
<crypko.objects.Crypko object at 0x7f1d39348e10>
>>> c.details.owner.username
'Bottersnike'
```

## Examples:

Examples can be found in `crypko/examples`.

