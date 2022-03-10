# Justin Biggs
#Lab 1 Assignment
USE company;

## 1) Retrieve first name and last name of all male employees with 
#salary more than 30000.

SELECT e.fname,e.lname
FROM employee AS e
WHERE sex = 'M' AND salary > 30000;

## 2) Retrieve locations of Research department projects.

SELECT DISTINCT p.plocation
FROM project AS p, department AS d, dept_locations AS dloc
WHERE d.dnumber=p.dnum AND d.dname = 'Research';

## 3) Retrieve first name, last name, and SSN of all 
#employees who work more than 9 hours on project #2.

SELECT DISTINCT e.fname,e.lname, e.ssn
FROM employee AS e, project AS p, works_on AS w
WHERE w.hours>9 AND w.essn=e.ssn AND p.pnumber=2 AND p.pnumber=w.pno;

## 4) Retrieve name, date of birth and relationship of all 
#female dependents of employees who work for department #5.

SELECT DISTINCT dep.dependent_name,dep.bdate,dep.relationship
FROM dependent AS dep, employee AS e
WHERE e.dno= 5 AND dep.sex = 'F' AND e.ssn=dep.essn;

## 5) Retrieve first name, last name and salary of employees who 
#manage departments with projects located in Houston.

SELECT DISTINCT e.fname, e.lname, e.salary
FROM employee AS e, project AS p, department AS d
WHERE p.plocation = 'Houston' AND p.dnum=d.dnumber AND d.mgr_ssn=e.ssn;