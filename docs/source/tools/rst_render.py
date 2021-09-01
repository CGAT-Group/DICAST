
def replace(text, replacewith, indent=0):
    return "{0}.. |{1}| replace:: {2}\n".format(" " * indent, text, replacewith)


def reference(text, target, indent=0):
    return "{0}.. _{1}: {2}\n".format(" " * indent, text, target)


def title(text, symbol="=", indent=0):
    underline = symbol * len(text)
    return "\n{0}{1}\n{0}{2}\n\n".format(" " * indent, text, underline)


def directive(directive, text, content, indent=0, addnewlines=True):
    content_out = indent_content(content, indent + 1, addnewlines)
    if addnewlines:
        sep = "\n"
    else:
        sep=""
    return "{0}.. {1}:: {2}{4}{0} {4}{3}{4}".format(" " * indent, directive, text, content_out, sep)


def note(content, indent=0):
    content_out = indent_content(content, indent + 1)
    return "{0}.. note::\n{0} \n{1}\n".format(" " * indent, content_out)

def warning(content, indent=0):
    content_out = indent_content(content, indent + 1)
    return "{0}.. warning::\n{0} \n{1}\n".format(" " * indent, content_out)

def indent_content(content, indent=1, addnewlines=True):
    content = content.split("\n")
    content_out = ""
    for line in content:
        content_out += "{0}{1}\n".format(" " * indent, line)
    if not addnewlines:
        content_out = content_out[:-2]
    return content_out


def codeblock(content, indent=0, language=None):
    content_out = indent_content(content, indent + 1)
    if language is None:
        return "{0}.. code-block::\n{0}\n{1}".format(" " * indent, content_out)
    else:
        return directive("code-block", language, content, indent)


def small_table(data, indent=0):
    colnum = len(data[0])
    colwidths = [0] * colnum
    for col in range(colnum):
        for row in data:
            if len(row[col]) > colwidths[col]:
                colwidths[col] = len(row[col])
    out = ""
    border = "  ".join(["=" * x for x in colwidths])
    out += "{0}{1}\n".format(" " * indent, border)
    for row in data:
        rowout = []
        for col in range(colnum):
            rowout.append("{0}{1}".format(row[col], " " * (colwidths[col] - len(row[col]))))
        out += "{0}{1}\n".format(" " * indent, "  ".join(rowout))
    out += "{0}{1}\n".format(" " * indent, border)
    return out


def bold(text, indent=0):
    return "{0}**{1}**".format(" " * indent, text)


def italic(text, indent=0):
    return "{0}*{1}*".format(" " * indent, text)


def comment(text, indent=0):
    return "{0}.. {1}\n".format(" " * indent, text)


def inline_link(text, link):
    return "`{0} <{1}>`_".format(text, link)

def doc_link(document):
    return f":doc:`{document}`"

def ref_link(text, link):
    return f":ref:`{text}<{link}>`"

def param_defin(parameter, content, indent=0):
    content_out = indent_content(content, indent + 1)
    return "{0}{1}\n{2}\n".format(" " * indent, parameter, content_out)


def li(content, bullet="*", indent=0):
    out = ""
    for line in content:
        out += "{0}{1} {2}\n".format(" " * indent, bullet, line)
    return out


def ul(content, bullet="#", indent=0):
    out = ""
    for line in content:
        out += "{0}{1}. {2}\n".format(" " * indent, bullet, line)
    return out
