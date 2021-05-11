<?xml version="1.0"?>
<!-- Was created with:
       xmlstarlet sel -C -t \
         -m "//*/Testing/Test[@Status='failed']" \
         -o '[ctest]  ' -v './/Name' -o ':' -n -n  \
         -v './/Results/Measurement/Value/text()' -n \
         -o '[ctest end]' -n -n \
         Testing/20210510-0943/Test.xml > ctest-error-output.xsl
-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exslt="http://exslt.org/common" version="1.0" extension-element-prefixes="exslt">
  <xsl:output omit-xml-declaration="yes" indent="no"/>
  <xsl:template match="/">
    <xsl:for-each select="//*/Testing/Test[@Status='failed']">
      <xsl:text>[ctest]  </xsl:text>
      <xsl:call-template name="value-of-template">
        <xsl:with-param name="select" select=".//Name"/>
      </xsl:call-template>
      <xsl:text>:</xsl:text>
      <xsl:value-of select="'&#10;'"/>
      <xsl:value-of select="'&#10;'"/>
      <xsl:call-template name="value-of-template">
        <xsl:with-param name="select" select=".//Results/Measurement/Value/text()"/>
      </xsl:call-template>
      <xsl:value-of select="'&#10;'"/>
      <xsl:text>[ctest end]</xsl:text>
      <xsl:value-of select="'&#10;'"/>
      <xsl:value-of select="'&#10;'"/>
    </xsl:for-each>
  </xsl:template>
  <xsl:template name="value-of-template">
    <xsl:param name="select"/>
    <xsl:value-of select="$select"/>
    <xsl:for-each select="exslt:node-set($select)[position()&gt;1]">
      <xsl:value-of select="'&#10;'"/>
      <xsl:value-of select="."/>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>
