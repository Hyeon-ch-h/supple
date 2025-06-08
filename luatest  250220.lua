-- FEMM script for extracting magnetic field data over a rectangular area


-- Value change section start
--
-- Define the rectangular area
x1 = 1.300    -- Starting x-coordinate
x2 = 2.700   -- Ending x-coordinate
y1 = 2.425    -- Starting y-coordinate
y2 = 3.825   -- Ending y-coordinate
dx = 0.001  -- Step size in x-direction
dy = 0.001  -- Step size in y-direction

--
-- Value change section end

-- Calculate the number of steps (assumes x2-x1 and y2-y1 are exactly divisible by dx, dy)
nx = (x2 - x1) / dx  -- for example: (10-0)/0.5 = 20
ny = (y2 - y1) / dy  -- for example: (10-0)/0.0002 = 50000

-- Open file for writing (FEMM environment)
filename = "magnetic_data.csv"
file = openfile(filename, "w")
if file == nil then
    print("Error: Unable to open file!")
    return
end

-- Write CSV header
write(file, "X,Y,Bx,By,Bmag\n")
-- Track last computed values
last_x = x1
last_y = y1

-- Use integer loop indices to compute x and y exactly
for i = 0, nx do
    local x_val = x1 + i * dx
    last_x = x_val  -- Update last processed x value

    for j = 0, ny do
        local y_val = y1 + j * dy
        last_y = y_val  -- Update last processed y value

        -- Get magnetic flux density components (Bx, By) using FEMM function
        A, Bx, By = mo_getpointvalues(x_val, y_val)

        -- Compute the magnitude of B: Bmag = sqrt(Bx^2 + By^2)
        Bmag = (Bx^2 + By^2) ^ 0.5

        -- Build a CSV-formatted line (concatenation automatically converts numbers to strings)
        local data_line = x_val .. "," .. y_val .. "," .. Bx .. "," .. By .. "," .. Bmag .. "\n"
        write(file, data_line)
    end

    if last_y ~= y2 then
        A, Bx, By = mo_getpointvalues(x_val, y2)
        Bmag = (Bx^2 + By^2) ^ 0.5
        local data_line = x_val .. "," .. y2 .. "," .. Bx .. "," .. By .. "," .. Bmag .. "\n"
        write(file, data_line)
    end

end

-- Ensure the last values (x2, y2) are included if they were skipped
if last_x ~= x2 then
    for j = 0, ny do
        local y_val = y1 + j * dy
        last_y = y_val  -- Update last processed y value
        A, Bx, By = mo_getpointvalues(x2, y_val)
        Bmag = (Bx^2 + By^2) ^ 0.5
        local data_line = x2 .. "," .. y_val .. "," .. Bx .. "," .. By .. "," .. Bmag .. "\n"
        write(file, data_line)
    end

    if last_y ~= y2 then
        A, Bx, By = mo_getpointvalues(x2, y2)
        Bmag = (Bx^2 + By^2) ^ 0.5
        local data_line = x2 .. "," .. y2 .. "," .. Bx .. "," .. By .. "," .. Bmag .. "\n"
        write(file, data_line)
    end

end

-- Close file
closefile(file)
print("Data saved successfully: " .. filename)