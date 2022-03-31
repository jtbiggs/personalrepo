#Make list of all project # for projects tht involve employee whose last name is smith either as a worker or a manager that controls the project.alter

USE company;

#smith is a worker
(SELECT w.pno
From works_on w, employee e
WHERE e.ssn = w.essn AND lname='Smith')
UNION
#Smith is manager
(SELECT p.pnumber
FROM project p, department d, employee e
WHERE p.dnum=d.dnumber AND d.mgr_ssn =e.ssn AND e.lname='Smith');

#Find all employees who address is in Houston
#Find all employees who was born in the 70s

SELECT fname,lname
FROM employee e
WHERE e.address LIKE '%Houston, TX';

SELECT fname,lname,bdate
FROM employee
WHERE bdate LIKE '__7_______';

#show resulting salaries if employees in project x are given a 10% raise

SELECT e.fname,e.lname,1.1*e.salary AS increased_sal
FROM employee e, works_on w, project p
WHERE e.ssn = w.essn AND w.pno=p.pnumber AND p.Pname='ProductX';

/* Query 6 - Retrieve a list of employees and the projects they are working on, ordered by department and, 
within each department, ordered alphabetically by last name, then first name.*/

SELECT DISTINCT d.dname,e.lname,e.fname,p.pname
FROM employee e, works_on w, project p, department d
WHERE e.ssn=w.essn AND w.pno=p.pnumber AND d.dnumber=e.dno
ORDER BY d.dname, e.lname, e.fname;
