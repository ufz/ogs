"""Script to test functions that are used to create the feature matrix."""

import pytest
from lxml import etree
from utils import apply_diff_to_base_xml


@pytest.fixture
def base_xml():
    return b"""<?xml version="1.0" encoding="UTF-8"?>
<root>
    <child1>value1</child1>
    <child2>value2</child2>
    <child3 age = "5"></child3>
    <child4 age = "6">value4</child4>
    <child5>
    <grandchild1 age = "4">value3</grandchild1>
    <grandchild2 age = "3">value4</grandchild2>
    </child5>
    <childn age = "10">value5</childn>
    <childn age = "9">value6</childn>
    <childn age = "8">value7</childn>
    <childn2 age = "11">value8</childn2>
    <childn2 age = "12">value9</childn2>
    <childn3 age = "13">value10</childn3>
    <childn3 age = "14">value11</childn3>
    <child6/>
</root>
"""


@pytest.fixture
def diff_xml_replace():
    return b"""<?xml version="1.0" encoding="ISO-8859-1"?>
<OpenGeoSysProjectDiff base_file="A2.prj">
    <replace sel="/*/child1"><child5>value10</child5></replace>
    <replace sel="/*/child2/text()">new_value2</replace>
    <replace sel = "/*/child3" type = "@age">7</replace>
    <replace sel = "/*/child4/@age">8</replace>
    <replace msel="/*/childn/text()">new_value3</replace>
    <replace sel="/*/childn2[2]/text()">new_value4</replace>
    <replace msel="/*/childn" type="@age">13</replace>
    <replace sel="/*/childn2[1]/@age">14</replace>
    <replace msel="/*/childn3"><childn4>new_value5</childn4></replace>
    <replace msel="/*/child5/*"><grandchildnew age = "2">new_value6</grandchildnew></replace>
    <replace sel = "/*/child6"/>
</OpenGeoSysProjectDiff>
"""


def test_process_diff_file_tags_replace(base_xml, diff_xml_replace):
    base_xml = etree.fromstring(base_xml).getroottree()
    diff_xml = etree.fromstring(diff_xml_replace).getroottree()
    base_xml = apply_diff_to_base_xml(diff_xml, base_xml, False)
    assert base_xml.xpath("/*/child1") == []
    assert base_xml.xpath("/*/child5/text()")[0] == "value10"
    assert base_xml.xpath("/*/child2/text()")[0] == "new_value2"
    assert base_xml.xpath("/*/child3/@age")[0] == "7"
    assert base_xml.xpath("/*/child4/@age")[0] == "8"

    assert base_xml.xpath("/*/childn[1]/text()")[0] == "new_value3"
    assert base_xml.xpath("/*/childn[2]/text()")[0] == "new_value3"
    assert base_xml.xpath("/*/childn[3]/text()")[0] == "new_value3"

    assert base_xml.xpath("/*/childn2[1]/text()")[0] == "value8"
    assert base_xml.xpath("/*/childn2[2]/text()")[0] == "new_value4"

    assert base_xml.xpath("/*/childn[1]/@age")[0] == "13"
    assert base_xml.xpath("/*/childn[2]/@age")[0] == "13"

    assert base_xml.xpath("/*/childn2[1]/@age")[0] == "14"
    assert base_xml.xpath("/*/childn2[2]/@age")[0] == "12"

    assert len(base_xml.xpath("/*/childn4")) == 2
    assert base_xml.xpath("/*/childn4[1]/text()")[0] == "new_value5"
    assert base_xml.xpath("/*/childn4[2]/text()")[0] == "new_value5"

    assert base_xml.xpath("/*/child5/grandchildnew[1]/@age")[0] == "2"
    assert base_xml.xpath("/*/child5/grandchildnew[2]/@age")[0] == "2"

    assert base_xml.xpath("/*/child5/grandchildnew[1]/text()")[0] == "new_value6"
    assert base_xml.xpath("/*/child5/grandchildnew[2]/text()")[0] == "new_value6"

    assert base_xml.xpath("/*/child6") == []


@pytest.fixture
def diff_xml_add():
    return b"""<?xml version="1.0" encoding="ISO-8859-1"?>
<OpenGeoSysProjectDiff base_file="A2.prj">
    <add sel="/*/child1" type = "@age">17</add>
    <add sel="/*/child2/@age">18</add>
    <add sel="/*/child3/text()">new_value5</add>
    <add sel="/*/child4/@age">8</add>
    <add msel="/*/childn/text()">new_value3</add>
    <add sel="/*/childn2[2]/text()">new_value4</add>
    <add msel="/*/childn" type="@age">13</add>
    <add sel="/*/childn2[1]/@age">14</add>
    <add msel="/*/childn3"><grandchild>new_value5</grandchild></add>
    <add msel="/*/child5/*"><greatgrandchild1 age="2">new_value6</greatgrandchild1><greatgrandchild2 age="1">new_value7</greatgrandchild2></add>
</OpenGeoSysProjectDiff>
"""


def test_process_diff_file_tags_add(base_xml, diff_xml_add):
    base_xml = etree.fromstring(base_xml).getroottree()
    diff_xml = etree.fromstring(diff_xml_add).getroottree()
    base_xml = apply_diff_to_base_xml(diff_xml, base_xml, False)

    # Test adding attribute (type=@age) - replaces existing attribute
    assert base_xml.xpath("/*/child3/text()")[0] == "new_value5"
    assert base_xml.xpath("/*/child1/@age")[0] == "17"
    assert base_xml.xpath("/*/child2/@age")[0] == "18"
    # Test adding to attribute (replaces existing attribute)
    assert base_xml.xpath("/*/child4/@age")[0] == "8"

    # Test multiple add on text (replaces existing text)
    assert base_xml.xpath("/*/childn[1]/text()")[0] == "new_value3"
    assert base_xml.xpath("/*/childn[2]/text()")[0] == "new_value3"
    assert base_xml.xpath("/*/childn[3]/text()")[0] == "new_value3"

    # Test single add on text (only second childn2 gets new value)
    assert base_xml.xpath("/*/childn2[1]/text()")[0] == "value8"
    assert base_xml.xpath("/*/childn2[2]/text()")[0] == "new_value4"

    # Test multiple add on attribute (replaces existing attribute)
    assert base_xml.xpath("/*/childn[1]/@age")[0] == "13"
    assert base_xml.xpath("/*/childn[2]/@age")[0] == "13"

    # Test single add on attribute (only first childn2 gets new value)
    assert base_xml.xpath("/*/childn2[1]/@age")[0] == "14"
    assert base_xml.xpath("/*/childn2[2]/@age")[0] == "12"

    # Test adding child elements (msel) - adds new children
    assert base_xml.xpath("/*/childn3[1]/grandchild/text()")[0] == "new_value5"
    assert base_xml.xpath("/*/childn3[2]/grandchild/text()")[0] == "new_value5"

    # Test adding multiple children (msel on child5/*)
    assert base_xml.xpath("/*/child5/grandchild1/greatgrandchild1/@age")[0] == "2"
    assert base_xml.xpath("/*/child5/grandchild2/greatgrandchild1/@age")[0] == "2"

    assert (
        base_xml.xpath("/*/child5/grandchild1/greatgrandchild1/text()")[0]
        == "new_value6"
    )
    assert (
        base_xml.xpath("/*/child5/grandchild2/greatgrandchild1/text()")[0]
        == "new_value6"
    )
    # Test adding multiple children (msel on child5/*)
    assert base_xml.xpath("/*/child5/grandchild1/greatgrandchild2/@age")[0] == "1"
    assert base_xml.xpath("/*/child5/grandchild2/greatgrandchild2/@age")[0] == "1"

    assert (
        base_xml.xpath("/*/child5/grandchild1/greatgrandchild2/text()")[0]
        == "new_value7"
    )
    assert (
        base_xml.xpath("/*/child5/grandchild2/greatgrandchild2/text()")[0]
        == "new_value7"
    )


@pytest.fixture
def diff_xml_remove():
    return b"""<?xml version="1.0" encoding="ISO-8859-1"?>
<OpenGeoSysProjectDiff base_file="A2.prj">
    <remove sel="/*/child1"/>
    <remove sel="/*/child2/text()"/>
    <remove sel = "/*/child3" type = "@age"/>
    <remove sel = "/*/child4/@age"/>
    <remove msel="/*/childn/text()"/>
    <remove sel="/*/childn2[2]/text()"/>
    <remove msel="/*/childn" type="@age"/>
    <remove sel="/*/childn2[1]/@age"/>
    <remove msel="/*/childn3"/>
</OpenGeoSysProjectDiff>
"""


def test_process_diff_file_tags_remove(base_xml, diff_xml_remove):
    base_xml = etree.fromstring(base_xml).getroottree()
    diff_xml = etree.fromstring(diff_xml_remove).getroottree()
    base_xml = apply_diff_to_base_xml(diff_xml, base_xml, False)

    # Test remove tag
    assert base_xml.xpath("/*/child1") == []

    # Test remove text content (set to empty string)
    assert base_xml.xpath("/*/child2/text()") == [""]

    # Test multiple remove on text
    assert base_xml.xpath("/*/childn[1]/text()") == [""]
    assert base_xml.xpath("/*/childn[2]/text()") == [""]
    assert base_xml.xpath("/*/childn[3]/text()") == [""]

    # Test remove attribute via xpath @age
    assert base_xml.xpath("/*/child4/@age") == []

    # Test single remove on text (only second childn2 gets value removed)
    assert base_xml.xpath("/*/childn2[1]/text()")[0] == "value8"
    assert base_xml.xpath("/*/childn2[2]/text()") == [""]

    # Test multiple remove on attribute (remove existing attribute)
    assert base_xml.xpath("/*/childn[1]/@age") == []
    assert base_xml.xpath("/*/childn[2]/@age") == []

    # Test single remove on attribute (only first childn2 loses its attribute)
    assert base_xml.xpath("/*/childn2[1]/@age") == []
    assert base_xml.xpath("/*/childn2[2]/@age")[0] == "12"


@pytest.fixture
def empty_selects():
    return [
        b"""<?xml version="1.0" encoding="ISO-8859-1"?>
    <OpenGeoSysProjectDiff base_file="A2.prj">"""
        + tag
        + b"""
    </OpenGeoSysProjectDiff>
    """
        for tag in [
            b'<replace sel="/*/chylde/@age">14</replace>',
            b'<replace sel="/*/chylde2"><chylde3></chylde3></replace>',
            b'<replace sel="/*/chylde4/text()">value</replace>',
            b'<filter sel="/*/child1"/>',
        ]
    ]


def test_empty_selector(base_xml, empty_selects):

    base_xml = etree.fromstring(base_xml).getroottree()

    # First 3 entries use valid-but-absent selectors → LookupError
    for empty_select in empty_selects[:3]:
        diff_xml = etree.fromstring(empty_select).getroottree()
        with pytest.raises(LookupError):
            base_xml = apply_diff_to_base_xml(diff_xml, base_xml, False)

    # 4th entry uses an unknown tag name → ValueError
    with pytest.raises(ValueError, match="Unknown tag: filter"):
        base_xml = apply_diff_to_base_xml(
            etree.fromstring(empty_selects[3]).getroottree(), base_xml, False
        )
