USE company;

# 1)
DELIMITER $$
CREATE TRIGGER new_dept
AFTER INSERT ON department
FOR EACH ROW 
BEGIN 
	INSERT INTO dept_locations
    VALUES(NEW.dnumber,'Houston');
END$$

/* TEST
INSERT INTO department
VALUES('HR', 2, 123456789, '2022-01-12');
DELETE FROM department
WHERE dname='HR';*/

# 2)
DELIMITER $$
CREATE TRIGGER salary_cap
BEFORE INSERT ON employee
FOR EACH ROW
BEGIN 
	DECLARE MSG VARCHAR(255);
    IF NEW.salary > (SELECT salary FROM employee e, department d WHERE NEW.dno=d.dnumber AND d.mgr_ssn=e.ssn)
    THEN 
    SET MSG = 'Supervisee salary is higher than supervisor salary.';
    SIGNAL SQLSTATE '45000' SET MESSAGE_TEXT = MSG;
    END IF;
END$$

/* TEST
DROP TRIGGER salary_cap;
INSERT INTO employee 
VALUES('Peter','', 'Parker', 123654987, '1996-03-11', 'Queens, New York City, New York', 'M', 30000, 888665555, 1);

DELETE FROM employee
WHERE ssn ='123654987';
*/

# 3) 
DELIMITER $$
CREATE TRIGGER new_employee_super
BEFORE INSERT ON employee
FOR EACH ROW
BEGIN 
	IF (NEW.super_ssn = '' OR NEW.super_ssn IS NULL)
    THEN
    SET NEW.super_ssn = (SELECT ssn FROM employee e, department d WHERE NEW.dno=d.dnumber AND d.mgr_ssn = e.ssn);
    END IF;
END$$

/* TEST 
INSERT INTO employee 
VALUES('Peter','', 'Parker', 123654987, '1996-03-11', 'Queens, New York City, New York', 'M', 30000, '', 5);
*/

# 4)
CREATE VIEW mgr_info
AS SELECT e.fname, e.lname, e.ssn, e.salary, e.dno, d.dname
FROM employee e, department d
WHERE e.ssn = d.mgr_ssn AND e.dno = d.dnumber; 

#5)

DROP VIEW project_view;
CREATE VIEW project_view(project_num, project_name, cont_dept_num, cont_dept_name, numbr_emp, total_proj_sal)
AS SELECT p.pnumber, p.pname, p.dnum, d.dname, COUNT(*), SUM(e.salary)
FROM employee e, project p, works_on w, department d
WHERE p.pnumber = w.pno AND w.essn=e.ssn AND p.dnum = d.dnumber
GROUP BY w.pno;