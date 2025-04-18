/*
 * Copyright (c) 2019, 2025, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.  Oracle designates this
 * particular file as subject to the "Classpath" exception as provided
 * by Oracle in the LICENSE file that accompanied this code.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 */

package jdk.jpackage.internal;

import java.text.MessageFormat;
import jdk.jpackage.internal.model.ConfigException;

final class PackageProperty {
    /**
     * Constructor
     *
     * @param name property name
     * @param expectedValue expected property value
     * @param substString substitution string to be placed in resource file to
     * be replaced with the expected property value by jpackage at package build
     * time
     * @param customResource name of custom resource from resource directory in
     * which this package property can be set
     */
    PackageProperty(String name, String expectedValue, String substString,
            String customResource) {
        this.name = name;
        this.expectedValue = expectedValue;
        this.substString = substString;
        this.customResource = customResource;
    }

    ConfigException verifyValue(String actualValue) {
        if (expectedValue.equals(actualValue)) {
            return null;
        }

        final String advice;
        if (substString != null) {
            advice = MessageFormat.format(I18N.getString(
                    "error.unexpected-package-property.advice"), substString,
                    actualValue, name, customResource);
        } else {
            advice = MessageFormat.format(I18N.getString(
                    "error.unexpected-default-package-property.advice"), name,
                    customResource);
        }

        return new ConfigException(MessageFormat.format(I18N.getString(
                "error.unexpected-package-property"), name,
                expectedValue, actualValue, customResource, substString), advice);
    }

    final String name;
    private final String expectedValue;
    private final String substString;
    private final String customResource;
}
