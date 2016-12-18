#!/bin/sh
rm $1.elki.txt
java -cp elki.jar de.lmu.ifi.dbs.elki.application.GeneratorXMLSpec -app.out $1.elki.txt -bymodel.spec $1.xml
