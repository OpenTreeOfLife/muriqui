import argparse, os, random, string
import muriqui as m
from StringIO import StringIO

MAX_ITEMS = 10
MAX_FLOAT = 100000.0
MAX_INT = 1000000000
MAX_KEY_LENGTH = 10
MAX_STRING_LENGTH = 10

simple_string_chars = string.letters + string.digits + string.punctuation

def get_random_string_utf8(length):

    # scour for utf-8 chars
    bytes = 8

    s = StringIO()
    i = 0
    while i < length:
        try:
            s.write(unicode(os.urandom(bytes), encoding='UTF-8'))
            i += 1
        except ValueError:
            continue

    return s.getvalue()

def get_random_string_ascii(length):
    s = StringIO()
    i = 0
    while i < length:
        s.write(random.sample(simple_string_chars,1)[0])
        i+=1
    return s.getvalue()

def get_random_float():
    return random.random() * random.randrange(MAX_FLOAT)

def get_random_int():
    return random.randrange(MAX_INT)

def get_random_value(depth):

    # halt deep recursion
    if depth > 8:
        return get_random_primitive()
        
    # generate containers infrequently
    r = random.randrange(8)
    if r == 0:
        value = []
        for i in range(random.randrange(MAX_ITEMS)):
            value.append(get_random_value())
    elif r == 1:
        value = {}
        for i in range(random.randrange(MAX_ITEMS)):
            value[get_random_string(MAX_KEY_LENGTH)] = get_random_value()
    else:
        # not a container, generate random primitive
        value = get_random_primitive()
    
    return value

def get_random_primitive():

    r = random.randrange(3)
    if r == 0:
        value = get_random_string(MAX_STRING_LENGTH)
    elif r == 1:
        value = get_random_float()
    else:
        value = get_random_int()

    return value

def get_random_body(number_of_elements):
    b = {}
    for i in range(number_of_elements):
        key = get_random_string(MAX_KEY_LENGTH)
        value = get_random_value()
        b[key] = value
    return b

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="generate random annotations")

    parser.add_argument("-l", "--label", nargs=1, required=True, help="The label to use for the output files.")
    
    parser.add_argument("-n", "--number-of-annotations", nargs=1, type=int, required=True, help="The number of random annotations to generate.")
    
    parser.add_argument("-s", "--random_number_seed", nargs=1, type=int, required=False, help="An integer random number seed. If none is supplied then the system time will be used.")

    parser.add_argument("-U", "--use-unicode", action="store_true", required=False, help="Whether to use unicode UTF-8 strings. If not specified, ascii strings will be generated.")
    
    args = parser.parse_args()
    label = args.label[0]
    count = args.number_of_annotations[0]
    if args.random_number_seed is not None:
        random.seed(args.random_number_seed[0])

    if args.use_unicode is not None and args.use_unicode == True:
        get_random_string = get_random_string_utf8 
    else:
        get_random_string = get_random_string_ascii

    # generate the annotations
    annotations = []
    for i in range(count):

        name = get_random_string(random.randrange(30))
        body = get_random_body(random.randrange(2,10))

        a = m.PhyloReferencedAnnotation()
        a.id = i
        a.annotated_at = "datetime"
        a.annotated_by = m.Entity.from_data({"name": name})
        a.body = body
                
        # set the target: generate random strings for ingroup
        a.target = m.ReferenceTarget.from_data({
            "type": "node",
            "included_ids": []
        })
        
        print a.summary
        exit()

#        node_leafset = set(n.leaf_nodes())
#        for l in tree.leaf_nodes():
#            if l not in node_leafset:
#                a.target.add_error_condition(m.CladeExcludesCheck(l.taxon.label))

        annotations.append(a)
    
    # write the tree to a file
    with open(label+".labeled.tre", "w") as treefile:
        treefile.write(tree.as_string("newick"))
    
    # write the annotations to a file
    with open(label+".annotations.json", "w") as datafile:        
        for a in annotations:
            datafile.write(a.summary + "\n")