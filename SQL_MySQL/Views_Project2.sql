#Create Views

USE mydb;
DROP VIEW relay_connections;
CREATE VIEW relay_connections(relay_SN, voltage_monitor_SN,breaker_number,current_monitor_SN)
AS SELECT r.`Serial Number`,v.`CCVT Serial Number`,b.`Breaker Number`,c.`CT Serial Number`
FROM `relay` AS r, `voltage monitor` v, breaker b, `current monitor` c
WHERE r.`serial number`=v.`relay serial number` AND v.`relay serial number`=r.`serial number` AND c.`relay serial number`=r.`serial number`;

#power transformer connections 
DROP VIEW transformer_connections;
CREATE VIEW transformer_connections
AS SELECT s.`name`, s.`ERCOT Code`, t.`serial number`, t.ratio,s.`power capacity`
FROM station s, (transformer t JOIN PT ON t.`serial number`=pt.`serial number`)
WHERE s.`ERCOT Code`=pt.`secondary side` AND s.`ERCOT Code`=t.`primary side`;

DROP VIEW number_of_breakers;
CREATE VIEW number_of_breakers
AS SELECT s.`name`, s.`ERCOT Code`,COUNT(*)
FROM station s, `relay` r, breaker b
WHERE s.`ERCOT Code`=r.`relay location` AND r.`serial number`=b.`relay serial number`
GROUP BY b.`relay serial number`;
