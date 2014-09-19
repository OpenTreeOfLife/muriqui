muriqui
=======

[![Build Status](https://secure.travis-ci.org/OpenTreeOfLife/muriqui.png)](http://travis-ci.org/OpenTreeOfLife/muriqui)

An annotation database for decorating nodes in trees

Couchdb application:
Couchdb views can be loaded using couchapp. 
---
sudo easy_install -U couchapp
---

Loading views

---
cd database/couchapp/ot
couchapp push . muriqui
---

In the above example, couchdb database is called ---muriqui---.

example query:
---
http://127.0.0.1:5984/muriqui/_design/couchapp/_view/by_date?startkey="2014/01/01"
---
