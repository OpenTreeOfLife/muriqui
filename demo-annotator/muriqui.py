import json

class GroupType:
    BRANCH, NODE, UNDEFINED = range(3)
    def to_str(c):
        if c == GroupType.BRANCH:
            return 'branch'
        elif c == GroupType.NODE:
            return 'node'
        assert c == GroupType.UNDEFINED
        return 'target type not defined'
    to_str = staticmethod(to_str)
    def to_code(c):
        if c.lower() == 'branch':
            return GroupType.BRANCH
        assert c.lower() == 'node'
        return GroupType.NODE
    to_code = staticmethod(to_code)

class PhyloReferencedAnnotation(object):
    def __init__(self):
        self.target = None
        self.annotated_at = None
        self.annotated_by = None
        self.body = None
        self.des = []
        self.exclude_ancs_of = []
        self.rooted_by = GroupType.UNDEFINED
        self.applied_to = []

    def set_annotated_by(self, entity):
        self.annotated_by = entity
    
    def set_annotated_at(self, datetime):
        self.annotated_at = datetime
    
    def set_target(self, reference):
        self.target = reference
        self.des = [str(i) for i in self.target['included_ids']]
        self.exclude_ancs_of = [str(i) for i in self.target.get('excluded_ids', [])]
        self.rooted_by = GroupType.to_code(self.target['type'])
    
    def set_body(self, body):
        self.body = body

    @classmethod
    def fromData(cls, data):
        a = cls()
        a.set_target(data['oa:hasTarget'])
        a.set_annotated_at(data['oa:annotatedAt'])
        a.set_annotated_by(data['oa:annotatedBy'])
        a.set_body(data['oa:hasBody'])
        a.applied_to = []
        return a

    def get_summary(self):
        return json.dumps(self.serialize())

    summary = property(get_summary)

    def serialize(self):
        return {'oa:hasTarget': self.target,
                'oa:annotatedBy': self.annotated_by,
                'oa:annotatedAt': self.annotated_at,
                'oa:hasBody': self.body,}
                