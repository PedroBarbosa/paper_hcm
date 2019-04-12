filters = [
    ('all', lambda x: x),
    ('coding', lambda x: x[ x['location'].str.match('coding')]),
    ('splicesite', lambda x: x[ x['location'].str.match('splicesite')]),
    ('5primeUTR', lambda x: x[ x['location'].str.match('5primeUTR')]),
    ('3primeUTR', lambda x: x[ x['location'].str.match('3primeUTR')]),
    ('deepintronic', lambda x: x[ x['location'].str.match('deepintronic')]),
    ('noncodingRNA', lambda x: x[ x['location'].str.match('noncodingRNA')]),
    ('mithocondrial', lambda x: x[x['location'].str.match('mithocondrial')]),
    ('unknown', lambda x: x[ x['location'].str.match('unknown')])
]
