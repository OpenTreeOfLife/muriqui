import argparse, datetime, os, random, string
import muriqui as m
from StringIO import StringIO

MAX_FLOAT = 100000.0
MAX_INT_VALUE = 1000000000

MAX_ITEMS = 10
MAX_KEY_LENGTH = 100
MAX_STRING_LENGTH = 10
MAX_TARGET_LENGTH = 100
MAX_ERROR_CONDITIONS = 10
MAX_WARNING_CONDITIONS = 10

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
    return random.randrange(MAX_INT_VALUE)

def get_random_value(depth=0):

    # halt deep recursion
    if depth > 8:
        return get_random_primitive()

    depth += 1
        
    # generate containers infrequently
    r = random.randrange(8)
    if r == 0:
        value = []
        for i in range(random.randrange(MAX_ITEMS)):
            value.append(get_random_value(depth))
    elif r == 1:
        value = {}
        for i in range(random.randrange(MAX_ITEMS)):
            value[get_random_string(MAX_KEY_LENGTH)] = get_random_value(depth)
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

def get_random_object(number_of_elements):
    b = {}
    for i in range(number_of_elements):
        key = get_random_string(MAX_KEY_LENGTH)
        value = get_random_value()
        b[key] = value
    return b
    
def get_random_condition(included_specifiers):
    i = random.randrange(2)
    # MonophylyCheck
    if i == 0: 
        specifiers = set()
        n = random.randrange(1,len(included_specifiers))
        c = m.MonophylyCheck(*random.sample(included_specifiers,n))
    # CladeExcludesCheck
    elif i == 1:
        specifiers = set()
        n = random.randrange(1,MAX_TARGET_LENGTH)
        while len(specifiers) < n:
            r = random.randrange(MAX_INT_VALUE)
            if r not in included_specifiers:
                specifiers.add(r)
        c = m.CladeExcludesCheck(*specifiers)
    return c

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

        e = m.Entity()
        e.name = get_random_string(random.randrange(30))
        e.url = "http://"+get_random_string(random.randrange(100))
        e.description = get_random_string(random.randrange(300))
        e.version = get_random_string(random.randrange(20))
        e.invocation = get_random_object(random.randrange(4))

        a = m.PhyloReferencedAnnotation()
        a.id = i
        a.annotated_at = datetime.datetime.now().isoformat()
        a.annotated_by = e
        a.body = get_random_object(random.randrange(2,10))
                
        # create a random reference type
        a.target = m.ReferenceTarget(random.sample(["node","branch"],1)[0])

        # generate random included ids        
        included = set()
        n_to_include = random.randrange(1,MAX_TARGET_LENGTH)
        while len(included) < n_to_include:
            included.add(random.randrange(MAX_INT_VALUE))

        # generate random excluded ids (only for branch references)
        if a.target.type == m.GroupType.BRANCH:
            excluded = set()
            n_to_exclude = random.randrange(MAX_TARGET_LENGTH)
            while len(excluded) < n_to_exclude:
                t = random.randrange(MAX_INT_VALUE)
                if t not in included:
                    excluded.add(t)
            a.target.exclude_specifiers(excluded)

        a.target.include_specifiers(included)

        # add zero or more error checks
        for i in range(random.randrange(MAX_ERROR_CONDITIONS)):
            a.target.add_error_condition(get_random_condition(included))
        
        # add zero or more warning checks
        for i in range(random.randrange(MAX_WARNING_CONDITIONS)):
            a.target.add_warning_condition(get_random_condition(included))
        
        print(a.summary)
        
        annotations.append(a)
        exit()
        
    # write the annotations to a file
    with open(label+".annotations.json", "w") as datafile:        
        for a in annotations:
            datafile.write(a.summary + "\n")