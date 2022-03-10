USE company;

#1) Retrieve first name and last name of employees whose birthday is in January.

SELECT fname,lname 
FROM employee AS e
WHERE bdate LIKE '_____01___';

/* 2) Show the resulting salaries if every employee working on the ‘ProductX’ 
project with a salary between $20000 and $40000 is given a 15% raise.*/

SELECT e.fname,e.lname,e.salary, 1.15*e.salary AS inc_salary
FROM employee AS e, project AS p, works_on AS w
WHERE p.pname = 'ProductX' AND p.pnumber=w.pno AND w.essn=e.ssn AND (e.salary BETWEEN 20000 AND 40000);

/* 3) List first name, last name and SSN of employees whose salary is less 
than the salary of any of the employees in department 4.*/

SELECT e.fname, e.lname, e.SSN
FROM employee e
WHERE e.salary < ANY (SELECT e.salary 
							FROM employee e
							WHERE dno = 4);

/* 4) Retrieve SSNs of all female employees who work on project numbers 10, 20, or 30.*/

SELECT DISTINCT e.SSN
FROM employee AS e, works_on AS w
WHERE e.sex = 'F' AND w.essn = e.ssn AND (w.essn IN (SELECT w.essn
								FROM WORKS_ON w WHERE pno=10)
				  OR  w.essn IN (SELECT w.essn
								FROM WORKS_ON w WHERE pno=20)
				  OR w.essn IN (SELECT w.essn
								FROM WORKS_ON w WHERE pno=30));

/* 5) For each project on which less than three employees work, retrieve project number, project name, 
and the average salary of employees who work on the project. */

SELECT p.pnumber, p.pname, AVG(e.salary) AS avg_salary
FROM project p, employee e, works_on w
WHERE p.pnumber = w.pno AND w.essn = e.ssn
GROUP BY p.pname, p.pnumber
HAVING COUNT(*) < 3;

								
