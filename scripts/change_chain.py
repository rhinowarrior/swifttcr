def change_chain(input_dir):
    #change I to D and J to E
    with open(input_dir, "r") as f:
        lines = f.readlines()
    with open(input_dir, "w") as f:
        for line in lines:
            if line.startswith("ATOM"):
                if line[21] == "I":
                    f.write(line.replace("I", "D"))
                elif line[21] == "J":
                    f.write(line.replace("J", "E"))
                else:
                    f.write(line)
            else:
                f.write(line)
