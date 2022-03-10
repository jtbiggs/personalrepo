INSERT INTO employee VALUES ('Richard','K','Marini','653298653','1962-12-30','98 Oak Forest, Kary, TX','M','37000','123456789',4);

DELETE FROM employee 
WHERE ssn = '653298653';

USE company;
DELETE FROM department
WHERE dname = 'HR';
INSERT INTO employee (fname,lname,ssn) VALUES ('Richard','Marini','653298653');

UPDATE employee
SET dno = 5
WHERE ssn = '653298653';

SELECT *
FROM employee
WHERE fname = 'John' AND minit='B' AND lname ='Smith';

#SELECT * FROM project
SELECT e.fname, e.lname, e.address FROM employee as e, department as d
WHERE e.dno = d.dnumber AND d.dname = 'Research';

SELECT DISTINCT fname, lname
FROM employee as e, department as d, project as p
WHERE d.mgr_ssn = e.ssn AND p.dnum=d.dnumber AND p.plocation = 'Stafford';