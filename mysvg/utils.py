def list_elements(svg_struct):
    finished = [];
    seen = [];
    seen.append(svg_struct);
    
    while (len(seen) > 0):
        elt = seen.pop();
        seen.extend(elt.__dict__.get('_subElements', []));
        finished.append(elt);

    return finished;

    
