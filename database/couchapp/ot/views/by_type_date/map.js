function(doc) {
  emit([doc["oa:hasBody"]["@type"],doc["oa:annotatedAt"]], doc);
}
