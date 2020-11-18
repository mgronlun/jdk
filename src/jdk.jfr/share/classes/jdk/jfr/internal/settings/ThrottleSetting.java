/*
 * Copyright (c) 2020, Oracle and/or its affiliates. All rights reserved.
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

package jdk.jfr.internal.settings;

import java.util.Objects;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import jdk.jfr.Description;
import jdk.jfr.Label;
import jdk.jfr.MetadataDefinition;
import jdk.jfr.Name;
import jdk.jfr.Timespan;
import jdk.jfr.internal.PlatformEventType;
import jdk.jfr.internal.Type;
import jdk.jfr.internal.Utils;

@MetadataDefinition
@Label("Event Emmission Throttle")
@Description("Event emissions throttle (events per second)")
@Name(Type.SETTINGS_PREFIX + "Throttle")
public final class ThrottleSetting extends JDKSettingControl {
    private final static long typeId = Type.getTypeId(ThrottleSetting.class);
    private static final String OFF = "off";

    private String value = "0";
    private final PlatformEventType eventType;

    public ThrottleSetting(PlatformEventType eventType) {
       this.eventType = Objects.requireNonNull(eventType);
    }

    @Override
    public String combine(Set<String> values) {
        long max = 0;
        String text = "0";
        for (String value : values) {
            long l = parseValueSafe(value);
            if (l == 0) {
                // throttling is disabled, accept everything
                return "0";
            }
            if (l > max) {
                text = value;
                max = l;
            }
        }
        return text;
    }

    @Override
    public void setValue(String value) {
        System.out.println("SetValue: " + value);
        long l =  parseValueSafe(value);
        this.value = value;
        eventType.setThrottle(l);
    }

    @Override
    public String getValue() {
        return value;
    }

    public static boolean isType(long typeId) {
        return ThrottleSetting.typeId == typeId;
    }

    public static long parseValueSafe(String value) {
        long result = -1L;
        try {
            result = Utils.parseThrottleValue(value);
        } catch (NumberFormatException nfe) {
            return -1L;
        }
        return result;
    }
    private static boolean isDisabled(String value) {
        return OFF.equals(value.trim().toLowerCase());
    }

    private static boolean isRate(String value) {
        // Expected input format is "x / y" or "x \ y" where x is a non-negative integer and
        // y is a time unit. Split the string at the delimiter.
        String[] result = value.split("[\\/\\\\]");
        return result.length == 2;
    }
    private static boolean isProbability(String value) {
        // Expected input format is "x %" or "0.x" where x is a non-negative integer and
        // y is a time unit. Split the string at the delimiter.
        String[] result = value.split("[\\/\\\\]");
        return result.length == 2;
    }
}
