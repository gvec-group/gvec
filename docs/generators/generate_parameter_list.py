import yaml
import textwrap
import re
from pathlib import Path


def generate_formatted_param(param, vals, formatting="markdown", enable_links=True):
    """
    generate a list of strings (each entry is a line)
    for a parameter named `param` and the dictionary `vals` which
    *  must contain 'description': a list of string
    *  can  contain 'subtitle'
    *  can contain 'required' True/False
    *  can  contain 'required_if': string
    *  must contain 'type': string
    *  must contain either 'allowed': string
    *               or 'allowed_table': list of (option, description) string tuple (each entry is a row in the table)
    *  can contain 'default': string

    formatting can be either "markdown" or "sceen", for now

    """

    # for linking between parameters (only for markdown formatting), the link inserting is done afterwards, using regexp.
    # simply add in the string expression
    #     "lnk_to_param(paramname)" , which means linktext=f"`{paramname}`
    # or  "lnk_to_param(linktext|paramname)"
    # and the whole "lnk_to_param(*)" is then replaced by the output of the call lnk_to_param(linktext,paramname)
    # if enable_links=False, or for the screen formatting, the whole "lnk_to_param(*)"" is replaced by only the `linktext`

    # helper functions for replacing occurences of "lnk_to_param(arg1)" or "lnk_to_param(arg1|arg2)"
    def lnk_to_param(link_text, param):
        if formatting == "markdown":
            if enable_links:
                return f"[{link_text}](#lnktoparam-{param})"
            else:
                return link_text
        if formatting == "screen":
            return link_text  # cannot enable links in screen output

    def replace_lnk_to_param(line):
        pattern = r"lnk_to_param\(([^)]+)\)"

        matches = re.findall(pattern, line)
        for match in matches:
            args = match.split("|")
            if len(args) == 1:
                repl = lnk_to_param(f"`{args[0]}`", args[0])
            if len(args) == 2:
                repl = lnk_to_param(args[0], args[1])
            line = line.replace(f"lnk_to_param({match})", repl)

        return line

    #####
    # main part
    #####

    outstr = []

    # debug, have open all dropdown open
    open_all = True

    if formatting == "markdown":
        outstr.append(f"::::{{dropdown}} `{param}`")
        if subtitle := vals.get("subtitle"):
            outstr[-1] += f" ({subtitle})"
        outstr.append(f":name: lnktoparam-{param}")
        if open_all:
            outstr.append(":open:")
        outstr.append("")
        for x in vals["description"]:
            outstr.append(replace_lnk_to_param(x))
            outstr.append("")

        if vals.get("required"):
            outstr.append("***required input!***")
            outstr.append("")

        if reqif := vals.get("required_if"):
            outstr.append(f"***required if:*** {replace_lnk_to_param(reqif)}")
            outstr.append("")

        outstr.append(f"**Type:** {vals['type']}")
        outstr.append("")

        if allowed := vals.get("allowed"):
            outstr.append(f"**Allowed Values:** {replace_lnk_to_param(allowed)}")
            outstr.append("")

        if table := vals.get("allowed_table"):
            outstr.append("**Allowed Values:**")
            # if any("■" in row[0] for row in table):
            #    outstr[-1] += " (■ default)"
            outstr.append(":::{list-table}")
            for row in table:
                outstr.append(f"*   - {replace_lnk_to_param(row[0])}")
                outstr.append(f"    - {replace_lnk_to_param(row[1])}")
            outstr.append(":::")  # list-table
            outstr.append("")

        if deflt := vals.get("default"):
            outstr.append(f"**Default:** {replace_lnk_to_param(deflt)}")
            outstr.append("")

        outstr.append("::::")  # dropdown

    if formatting == "screen":
        # helper for text wrap:
        def wrap_swid(pre, line, linelen):
            maxlen = linelen - len(pre)
            if len(line) <= maxlen:
                return [pre + line]
            else:
                wline = textwrap.wrap(line, width=maxlen)
                out = []
                out.append(pre + wline[0])
                [out.append(" " * len(pre) + wl) for wl in wline[1:]]
                return out

        swid = 80  # screen witdh

        outstr.append("┌" + "─" * (swid + 2) + "┐")
        sepline = "├" + "─" * (swid + 2) + "┤"

        line = f"`{param}`"
        if subtitle := vals.get("subtitle"):
            line += f" ({subtitle})"
        outstr.append(f"| {line:{swid}} |")
        outstr.append(sepline)

        descr = [replace_lnk_to_param(x) for x in vals["description"]]
        wdescr = []
        for x in descr:
            wdescr += wrap_swid("   ", x, swid)
        [outstr.append(f"| {x:{swid}} |") for x in wdescr]

        if vals.get("required"):
            outstr.append("├" + "─" * (swid + 2) + "┤")
            line = "required input!"
            outstr.append(f"|    {line:{swid-3}} |")
        if reqif := vals.get("required_if"):
            outstr.append("├" + "─" * (swid + 2) + "┤")
            wreqif = wrap_swid("   required if: ", replace_lnk_to_param(reqif), swid)
            [outstr.append(f"| {x:{swid}} |") for x in wreqif]
            # reqif_swid = swid-3-13
            # if len(reqif_repl) <= reqif_swid:
            #    outstr.append(f"|    required if: {reqif_repl:{reqif_swid}} |")
            # else:
            #    w_reqif = textwrap.wrap(reqif_repl, width=reqif_swid)
            #    outstr.append(f"|    required if: {w_reqif[0]:{reqif_swid}} |")
            #    for wr in w_reqif[1:]:
            #        outstr.append(f"|                 {wr:{reqif_swid}} |")

        outstr.append(sepline)
        line = f"Type: {vals['type']}"
        outstr.append(f"|    {line:{swid-3}} |")

        if alwd := vals.get("allowed"):
            outstr.append(sepline)
            walwd = wrap_swid("   Allowed Values: ", replace_lnk_to_param(alwd), swid)
            [outstr.append(f"| {x:{swid}} |") for x in walwd]

        if table := vals.get("allowed_table"):
            outstr.append(sepline)
            line = "Allowed Values:"
            outstr.append(f"|    {line:{swid-3}} |")

            maxcol1 = max([len(row[0]) for row in table])
            maxcol2_swid = swid - 10 - maxcol1

            nrows = len(table)
            for irow, row in enumerate(table):
                if irow == 0:
                    outstr.append(
                        (
                            "|    "
                            + "┌"
                            + "─" * (maxcol1 + 2)
                            + "┬"
                            + "─" * (maxcol2_swid + 2)
                            + "┐ |"
                        )
                    )
                col1 = replace_lnk_to_param(row[0])
                wcol2 = wrap_swid("", replace_lnk_to_param(row[1]), maxcol2_swid)
                for jrow, wc in enumerate(wcol2):
                    if jrow == 0:
                        outstr.append(
                            "|    " + f"| {col1:{maxcol1}} | {wc:{maxcol2_swid}} | |"
                        )
                    else:
                        outstr.append(
                            "|    | " + " " * maxcol1 + f" | {wc:{maxcol2_swid}} | |"
                        )
                if irow < nrows - 1:
                    outstr.append(
                        "|    "
                        + "├"
                        + "─" * (maxcol1 + 2)
                        + "┼"
                        + "─" * (maxcol2_swid + 2)
                        + "┤ |"
                    )
                if irow == nrows - 1:
                    outstr.append(
                        "|    "
                        + "└"
                        + "─" * (maxcol1 + 2)
                        + "┴"
                        + "─" * (maxcol2_swid + 2)
                        + "┘ |"
                    )

        if dflt := vals.get("default"):
            outstr.append(sepline)
            wdflt = wrap_swid("   Default: ", replace_lnk_to_param(dflt), swid)
            [outstr.append(f"| {x:{swid}} |") for x in wdflt]

        outstr.append("└" + "─" * (swid + 2) + "┘")
    outstr.append("")

    return outstr


def format_parameter_list(
    yamlfile, formatting="markdown", output_file=None, filter={}, enable_links=True
):
    """
    filter is a dictionary, with "name":searchstring and/or "category":searchstring  as possibilities.
    This then only shows the entrys where the searchstring is found (ignore case) in the parameter name  and/or in the category  of the parameter
    """
    with open(yamlfile, "r") as f:
        dict = yaml.safe_load(f)
    out = []

    for param, vals in dict.items():
        if filt_name := filter.get("name"):
            if filt_name.lower() not in param.lower():
                continue
        if filt_category := filter.get("category"):
            if not any(
                [(filt_category.lower() in x.lower()) for x in vals["category"]]
            ):
                continue
        out += generate_formatted_param(
            param, vals, formatting=formatting, enable_links=enable_links
        )

    if output_file:
        with open(output_file, "w") as f:
            for line in out:
                f.write(line + "\n")
    else:
        for line in out:
            print(line)
