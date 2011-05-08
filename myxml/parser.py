from lxml import etree

nsmap = {
        'sodipodi': 'http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd',
        'cc': 'http://web.resource.org/cc/',
        'svg': 'http://www.w3.org/2000/svg',
        'dc': 'http://purl.org/dc/elements/1.1/',
        'xlink': 'http://www.w3.org/1999/xlink',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'inkscape': 'http://www.inkscape.org/namespaces/inkscape'
        }

    
def parse(data):
    assert type(data) == str, 'Data type should be string'
    data = etree.XML(data)
    return data

def svg_elt(tree, dtype):
    return tree.xpath('//svg:{0}'.format(dtype),namespaces=nsmap)

'''
# All svg text elements
>>> data.xpath('//svg:text',namespaces=nsmap)
[<Element {http://www.w3.org/2000/svg}text at b7cfc9dc>]
# All svg text elements with id="libcode-00"
>>> data.xpath('//svg:text[@id="libcode-00"]',namespaces=nsmap)
[<Element {http://www.w3.org/2000/svg}text at b7cfc9dc>]
# TSPAN child elements of text elements with id="libcode-00"
>>> data.xpath('//svg:text[@id="libcode-00"]/svg:tspan',namespaces=nsmap)
[<Element {http://www.w3.org/2000/svg}tspan at b7cfc964>]
# All text elements with id starting with "libcode"
>>> data.xpath('//svg:text[fn:startswith(@id,"libcode")]',namespaces=nsmap)
[<Element {http://www.w3.org/2000/svg}text at b7cfcc34>]
# Iterate text elements, access tspan child
>>> for elem in data.xpath('//svg:text[fn:startswith(@id,"libcode")]',namespaces=nsmap):
...     tp = elem.xpath('./svg:tspan',namespaces=nsmap)[0]
...     tp.text = "new text"
'''
